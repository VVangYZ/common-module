import numpy as np
import pandas as pd

E = 2.06 *10 ** 5
type_1 = [
    [1, 1],
    [[1, 2],[-1, -2]],
    [-1, -2], 
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 3],
    [3, 3]
]

type_2 = [
    [[2, 3], [3, 4]],
    [[2, 2], [3, 4]],
    [[2, 2], [3, 3]]
]

alpha = [
    [0.41, 0.986, 0.152],
    [0.65, 0.965, 0.3],
    [[0.73, 0.906, 0.595], [0.73, 1.216, 0.302]],
    [[1.35, 0.868, 0.915], [1.35, 1.375, 0.432]]
]

class Steel_c():
    """受压构件验算

    参数：
    l：杆件长度（单位：m）
    i：回转半径（单位：cm）
    u：计算长度系数
    fy：屈服强度（单位：MPa）
    """

    def __init__(self, l, i, u, fy=345):
        self.i = np.array(i)
        self.u = np.array(u)
        self.fy = fy
        self.E = E
        self.l = l * self.u * 100
        self.lam = self.l / self.i
        self.lam_b = self.lam / np.pi * np.sqrt(self.fy / self.E)
    
    def get_phi(self, t, shape, b=1, h=1, n=0):
        """获取截面类型，得到相关系数

        参数：
        t：板厚（单位：mm）
        shape：杆件截面形状，取 P85 页表格列数
        b：工字型截面宽度（单位：mm）
        h：工字型截面高度（单位：mm）
        n：如果为第二个表格，取第几小行
        """

        # 获取截面类型
        self.t = t
        if t < 40:
            if shape == 2:
                self.type = np.array(type_1[1][int(b / h > 0.8)])
            else:
                self.type = np.array(type_1[shape - 1])
        else:
            self.type = np.array(type_2[shape - 1][n])
        
        if self.type[0] < 0:
            if self.fy > 235:
                self.type = - self.type
            else:
                self.type = - self.type + 1

        # 获取系数
        self.alpha = []
        for i, j in enumerate(self.type):
            if j <= 2:
                self.alpha.append(alpha[j - 1])
            else:
                m = int(self.lam_b[i] > 1.05)
                self.alpha.append(alpha[j - 1][m])
        
        self.alpha = np.array(self.alpha)
        
        # 计算轴心压杆稳定系数
        phi_1 = 1 - self.alpha[:, 0] * self.lam_b ** 2
        phi_2 = 1 / (2 * self.lam_b ** 2) * \
            ((self.alpha[:, 1] + self.alpha[:, 2] * self.lam_b + self.lam_b ** 2)\
             - ((self.alpha[:, 1] + self.alpha[:, 2] * self.lam_b + self.lam_b ** 2)\
             ** 2 - 4 * self.lam_b ** 2) ** 0.5)
        
        self.phi = np.where(self.lam_b > 0.215, phi_2, phi_1)

        return self.phi


    def get_ncrd(self, A, fd=305):
        self.A = A
        self.fd = fd
        self.Ncrd = min(self.phi) * A * fd * 100 / 1000

        return self.Ncrd


        

        

