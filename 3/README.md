**利用了线性代数Eigen库作为辅助完成矩阵乘法等基本操作：**
- 其中Matrix<double, Dynamic, Dynamic>声明了一个元素为double的动态矩阵，其大小可运行时指定。
- MatrixXd为double类型的动态矩阵
- MatrixXcd为复数类型的动态矩阵
- 对于矩阵A，可直接用A(i,j)来访问第i行第j列的元素
- VectorXd为动态向量，其长度可变