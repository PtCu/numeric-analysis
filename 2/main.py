import csv
import numpy as np
from numpy.core.fromnumeric import transpose
from numpy.core.multiarray import concatenate


#男生身高
m_H = np.loadtxt(open("numeric_analytics\lab2\man.csv",
                      encoding="utf-8"), delimiter=",", skiprows=1, usecols=(1))
#男生体重
m_W = np.loadtxt(open("numeric_analytics\lab2\man.csv",
                      encoding="utf-8"), delimiter=",", skiprows=1, usecols=(2))
#女生身高                      
w_H = np.loadtxt(open("numeric_analytics\lab2\woman.csv",
                      encoding="utf-8"), delimiter=",", skiprows=1, usecols=(1))
#女生体重
w_W = np.loadtxt(open("numeric_analytics\lab2\woman.csv",
                      encoding="utf-8"), delimiter=",", skiprows=1, usecols=(2))
#总身高
t_H = np.concatenate((m_H, w_H), axis=0)
#总体重
t_W = np.concatenate((m_W, w_W), axis=0)

# gram-schmidt正交化
def gram_schmidt(A):
    Q = np.zeros_like(A)
    j = 0
    for a in A.T:  # 对每一个向量a
        y = np.copy(a)
        for i in range(0, j):  # 对于每个其他的维度
            y -= np.dot(np.dot(Q[:, i].T, a), Q[:, i])  # y减去在a上投影，得到正交的向量
        e = y / np.linalg.norm(y)  # 归一化
        Q[:, j] = e
        j += 1
    R = np.dot(Q.T, A)
    return (Q, R)


# 求最小二乘解
def LS(_A, _b):
    b = np.matmul(transpose(_A), _b)
    # 求出QR分解
    (Q, R) = gram_schmidt(np.matmul(transpose(_A), _A))
    # 变换等式求最小二乘解
    x = np.linalg.solve(R, np.matmul(Q.T, b))
    return x

#解甲、乙、丁
def solveNormal(_power, H, W):
    # 求平均值、方差、中位数
    C = W[:]/np.power(H[:], _power)

    print("平均值为："+str(np.mean(C)))
    print("方差为："+str(np.var(C)))
    print("中位数为："+str(np.median(C)))

    # 求最小二乘解
    A = np.zeros((len(H), 1))
    A[:, 0] = np.power(H[:], _power)
    b = np.zeros((len(W), 1))
    b[:, 0] = W[:]
    x = LS(A, b)

    c = x[0][0]
    print("最小二乘解为: "+str(c))

#解丙
def solveAbNormal(H, W):
    A = np.ones((len(H), 2))
    A[:, 0] = 1
    A[:, 1] = np.log(H[:])
    b = np.zeros((len(W), 1))
    b[:, 0] = np.log(W[:])

    x = LS(A, b)
    k = x[0][0]
    c2 = x[1][0]
    c1 = np.exp(k)
    print("情况丙：")
    print("c1为；"+str(c1)+"  c2为："+str(c2))
    return c1, c2

# 求解情况戊
def solveFive(c1, H, W):
     # 求平均值、方差、中位数
    C = (np.log(W[:])-np.log(c1))/np.log(H[:])

    print("平均值为："+str(np.mean(C)))
    print("方差为："+str(np.var(C)))
    print("中位数为："+str(np.median(C)))

    b = np.zeros(len(W), 1)
    b[:, 0] = np.log(W[:])-np.log(c1)
    A = np.zeros(len(H), 1)
    A[:, 0] = np.log(H[:])

    c = LS(A, b)
    print("最小二乘解为: "+str(c))


def solve():
    print("甲\n")
    print("男生")
    solveNormal(2, m_H, m_W)
    print("女生")
    solveNormal(2, w_H, w_W)
    print("总体")
    solveNormal(2, t_H, t_W)

    print("\n")
    print("乙\n")
    print("男生")
    solveNormal(3, m_H, m_W)
    print("女生")
    solveNormal(3, w_H, w_W)
    print("总体")
    solveNormal(3, t_H, t_W)

    print("\n")
    print("丙\n")
    print("男生")
    m_c1, m_c2 = solveAbNormal(m_H, m_W)
    print("女生")
    w_c1, w_c2 = solveAbNormal(w_H, w_W)
    print("总体")
    t_c1, t_c2 = solveAbNormal(t_H, t_W)

    print("\n")
    print("丁\n")
    print("男生")
    solveNormal(m_c2, m_H, m_W)
    print("女生")
    solveNormal(w_c2, w_H, w_W)
    print("总体")
    solveNormal(t_c2, t_H, t_W)

    print("\n")
    print("戊\n")
    print("男生")
    solveNormal(m_c1, m_H, m_W)
    print("女生")
    solveNormal(w_c1, w_H, w_W)
    print("总体")
    solveNormal(t_c1, t_H, t_W)


if __name__ == "__main__":
    solve()
