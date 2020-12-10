#include <iostream>
#include "math.hpp"

int main()
{

    Matrix<double, 8, 8> A;
    Matrix<double, 10, 10> B;
    Matrix<double, 12, 12> C;
    A << 611, 196, -192, 407, -8, -52, -49, 29,
        196, 899, 113, -192, -71, -43, -8, -44,
        -192, 113, 899, 196, 61, 49, 8, 52,
        407, -192, 196, 611, 8, 44, 59, -23,
        -8, -71, 61, 8, 411, -599, 208, 208,
        -52, -43, 49, 44, -599, 411, 208, 208,
        -49, -8, 8, 59, 208, 208, 99, -911,
        29, -44, 52, -23, 208, 208, -911, 99,

        B << 1/1,  1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8,  1/9,  1/10,
	1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8,  1/9,  1/10, 1/11,
	1/3,  1/4,  1/5,  1/6,  1/7,  1/8,  1/9,  1/10, 1/11, 1/12,
	1/4,  1/5,  1/6,  1/7,  1/8,  1/9,  1/10, 1/11, 1/12, 1/13,
	1/5,  1/6,  1/7,  1/8,  1/9,  1/10, 1/11, 1/12, 1/13, 1/14,
	1/6,  1/7,  1/8,  1/9,  1/10, 1/11, 1/12, 1/13, 1/14, 1/15,
	1/7,  1/8,  1/9,  1/10, 1/11, 1/12, 1/13, 1/14, 1/15, 1/16,
	1/8,  1/9,  1/10, 1/11, 1/12, 1/13, 1/14, 1/15, 1/16, 1/17,
	1/9,  1/10, 1/11, 1/12, 1/13, 1/14, 1/15, 1/16, 1/17, 1/18,
	1/10, 1/11, 1/12, 1/13, 1/14, 1/15, 1/16, 1/17, 1/18, 1/19,

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
            D(i, j) = sqrt(2.0f / 21.0f) * sin((double)i * j * acos(-1.0) / 21);
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

    /*A*/
    //幂法
    double pld;
    std::vector<double> env;
    power_eng(A, pld, env);
    std::cout<<"幂法"<<std::endl;
    std::cout << "A的特征向量" << std::endl;
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl<< "A的特征值   " <<pld << std::endl;

    env.clear();
    jacobi_eng(A,env);
    std::cout<<"jacobi"<<std::endl;
    std::cout << "A的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl;
    env.clear();


    unshifted_qr(A,env);
    std::cout<<"无移动QR"<<std::endl;
    std::cout << "A的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
 std::cout<<std::endl;
 env.clear();

    std::vector<std::complex<double>> env1;
     shifted_qr(A,env1);
     std::cout<<"移动QR"<<std::endl;
        std::cout << "A的特征值" << std::endl;
    sort(env1.begin(),env1.end(),[](const std::complex<double>& a,const std::complex<double>& b){return a.real()<b.real();});
    for (size_t i = 0; i < env1.size(); ++i)
    {
        std::cout << env1[i] <<" ";
    }
    std::cout<<std::endl;
      env1.clear();
   
    /*B*/

    power_eng(B, pld, env);
    std::cout<<"幂法"<<std::endl;
    std::cout << "B的特征向量" << std::endl;
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl<< "B的特征值   " <<pld << std::endl;
    env.clear();

    jacobi_eng(B,env);
    std::cout<<"jacobi"<<std::endl;
    std::cout << "B的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl;
    env.clear();


    unshifted_qr(B,env);
    std::cout<<"无移动QR"<<std::endl;
    std::cout << "B的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
 std::cout<<std::endl;
env.clear();
   
     shifted_qr(A,env1);
     std::cout<<"移动QR"<<std::endl;
        std::cout << "B的特征值" << std::endl;
    sort(env1.begin(),env1.end(),[](const std::complex<double>& a,const std::complex<double>& b){return a.real()<b.real();});
    for (size_t i = 0; i < env1.size(); ++i)
    {
        std::cout << env1[i] <<" ";
    }
    std::cout<<std::endl;
     env1.clear();



    //c
    power_eng(C, pld, env);
    std::cout<<"幂法"<<std::endl;
    std::cout << "C的特征向量" << std::endl;
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl<< "C的特征值   " <<pld << std::endl;
    env.clear();
  
    jacobi_eng(C,env);
    std::cout<<"jacobi"<<std::endl;
    std::cout << "C的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl;
    env.clear();


    unshifted_qr(C,env);
    std::cout<<"无移动QR"<<std::endl;
    std::cout << "C的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
 std::cout<<std::endl;
  env.clear();
    
     shifted_qr(C,env1);
     std::cout<<"移动QR"<<std::endl;
        std::cout << "C的特征值" << std::endl;
    sort(env1.begin(),env1.end(),[](const std::complex<double>& a,const std::complex<double>& b){return a.real()<b.real();});
    for (size_t i = 0; i < env1.size(); ++i)
    {
        std::cout << env1[i] <<" ";
    }
    std::cout<<std::endl;
      env1.clear();


        //d
    power_eng(D, pld, env);
    std::cout<<"幂法"<<std::endl;
    std::cout << "D的特征向量" << std::endl;
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl<< "D的特征值   " <<pld << std::endl;
 env.clear();
   
    jacobi_eng(D,env);
    std::cout<<"jacobi"<<std::endl;
    std::cout << "D的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl;
    env.clear();


    unshifted_qr(D,env);
    std::cout<<"无移动QR"<<std::endl;
    std::cout << "D的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
 std::cout<<std::endl;
 env.clear();
   
     shifted_qr(D,env1);
     std::cout<<"移动QR"<<std::endl;
        std::cout << "C的特征值" << std::endl;
    sort(env1.begin(),env1.end(),[](const std::complex<double>& a,const std::complex<double>& b){return a.real()<b.real();});
    for (size_t i = 0; i < env1.size(); ++i)
    {
        std::cout << env1[i] <<" ";
    }
    std::cout<<std::endl;
    env1.clear();

    //e
     power_eng(E, pld, env);
    std::cout<<"幂法"<<std::endl;
    std::cout << "E的特征向量" << std::endl;
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl<< "E的特征值   " <<pld << std::endl;
 env.clear();
 
    jacobi_eng(E,env);
    std::cout<<"jacobi"<<std::endl;
    std::cout << "E的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
    std::cout << std::endl;
    env.clear();


    unshifted_qr(E,env);
    std::cout<<"无移动QR"<<std::endl;
    std::cout << "E的特征值" << std::endl;
    sort(env.begin(),env.end());
    for (size_t i = 0; i < env.size(); ++i)
    {
        std::cout << env[i] <<" ";
    }
 std::cout<<std::endl;
   env.clear();
    
     shifted_qr(E,env1);
     std::cout<<"移动QR"<<std::endl;
        std::cout << "E的特征值" << std::endl;
    sort(env1.begin(),env1.end(),[](const std::complex<double>& a,const std::complex<double>& b){return a.real()<b.real();});
    for (size_t i = 0; i < env1.size(); ++i)
    {
        std::cout << env1[i] <<" ";
    }
    std::cout<<std::endl;
    env1.clear();
    return 0;
}
