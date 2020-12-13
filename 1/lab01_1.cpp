#include <bits/stdc++.h>
using namespace std;

//单精度浮点数求平方根算法
float q_sqrt(float c)
{
    float c_half = 0.5f * c;
    int i = *(int *)&c; //获取float类型的位数
    i = 0x5f375a86 - (i >> 1);
    c = *(float *)&i;
    c = c * (1.5f - c_half * c * c);
    c = c * (1.5f - c_half * c * c);
    return 1 / c;
}

//将一个数开跟并取倒数
float invSqrt(float x)
{
    float xhalf = 0.5 * x;
    int i = *(int *)&x;            // get bits for floating value
    i = 0x5f3759df - (i >> 1);     // gives initial guess
    x = *(float *)&i;              // convert bits back to float
    x = x * (1.5 - xhalf * x * x); // Newton step
    return x;
}
int main(){
    q_sqrt(4.0);
}