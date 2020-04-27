"""用于进行桥梁静阵风荷载的计算脚本

基于《公路桥梁抗风设计规范》（JTG/T 3360-01-2018）前四章进行计算
"""

import pandas as pd
import numpy as np
import os

current_path = os.path.dirname(__file__)


# 一些常数
l_alpha0 = {'A': 0.12, 'B': 0.16, 'C': 0.22, 'D': 0.30}
l_z0 = {'A': 0.01, 'B': 0.05, 'C': 0.3, 'D': 1.0}
l_kc = {'A': 1.174, 'B': 1.0, 'C': 0.785, 'D': 0.564}
l_kf = [1.05, 1.02, 1.00]
l_ksf = [[0.88, 0.84, 0.78], [0.92, 0.88, 0.84]]
l_GV_l = pd.read_csv(current_path + r'\wind\GV_l.csv')
l_GV_h = pd.read_csv(current_path + r'\wind\GV_h.csv')
l_CD = pd.read_csv(current_path + r'\wind\CD.csv')
l_CH = pd.read_csv(current_path + r'\wind\CH.csv')
l_eta = pd.read_csv(current_path + r'\wind\eta.csv')


class Wind():
    """定义桥位整体风效应情况
    
    参数：
        r：危险区域(1-3)
        sc：地表类别(A-D)
        w1：十年风速(m/s)
        w2：百年风速(m/s)
        year：施工年限
    """

    def __init__(self, r=2, sc='B', w1=0, w2=0, year=1):
        """构造函数"""
        self.R = r
        self.surface_class = sc
        self.W1 = min(w1, 25)
        self.W2 = w2
        self.U10 = max(24.5, w2)

        self.alpha0 = l_alpha0[sc]  # 地表粗糙度系数（4.2.1）
        self.z0 = l_z0[sc]          # 地表粗糙高度（4.2.1）
        self.kc = l_kc[sc]          # 地表类别转换系数 (4.2.4)
        self.Us10 = self.kc * self.U10
        self.kf = l_kf[r-1]         # 抗风风险系数（4.2.6）
        self.ksf = l_ksf[year<=3][r-1]  # 施工抗风风险系数（4.2.9）
    
    def get_situation(self):
        """区域地表参数"""
        print(f'桥梁抗风风险区域：R{self.R}')
        print(f'地表分类：{self.surface_class}')

    def get_U(self):
        """各风速计算结果"""
        print(f'十年重现期风作用水平 W1：{self.W1:.3f}')
        print(f'百年重现期风作用水平 W2：{self.W2:.3f}')
        print(f'基本风速 U10：{self.U10:.3f} m/s')
        print(f'桥梁设计基准风速为 Us10：{self.Us10:.3f} m/s')
    
    def get_k(self):
        """各参数"""
        print(f'地表粗糙度系数（4.2.1）：{self.alpha0:.3f}')
        print(f'地表粗糙高度（4.2.1）：{self.z0:.3f} m')
        print(f'地表类别转换系数 (4.2.4)：{self.kc:.3f}')
        print(f'抗风风险系数（4.2.6）：{self.kf:.3f}')
        print(f'施工抗风风险系数（4.2.9）:{self.ksf:.3f}')


class BeamWind():
    def __init__(self, wind, shape, z, l, rho=1.25, ud=25):
        """主梁风荷载
            
            参数：
                wind：设计风实例，如缺少则输入 0，直接输入 ud
                shape：主梁断面形状（list，单位：m）
                    -闭口流线型箱梁：[高度]
                    -工字形梁和Ⅱ字梁：[宽度，高度]
                    -箱型：[宽度，高度，腹板倾角]
                    -桁架：[杆件形状，实面积比，中心距，高度，*各杆件高度]
                        --杆件形状：1、2（矩形或 H 形、圆形）
                z：主梁高度(m)
                l：单跨长度(m)
                rho：空气密度(kg/m3)
                ud：如果缺少设计风资料，用于估计风速值，默认地表为 A (m/s)
            """
        self.wind = wind
        self.rho = rho
        self.z = z
        self.l = l
        self.shape = len(shape)
        if wind == 0:
            self.Ud = ud
            self.GV = np.interp(l, l_GV_l['L'], l_GV_l['A'])
        else:
            self.Ud = wind.kf * (z / 10) ** wind.alpha0 * wind.Us10
            self.GV = np.interp(l, l_GV_l['L'], l_GV_l[wind.surface_class])
        
        self.Ug = self.GV * self.Ud

        if len(shape) == 1:
            # 闭口流线型箱梁，参数为高度
            self.D = shape[0]
            self.CH = 1.1
            self.Fg = 1/2 * rho * self.Ug ** 2 * self.CH * z
        elif len(shape) == 2:
            # 工字形或Ⅱ字形，参数为宽度与高度
            self.D = shape[1]
            bd = shape[0] / shape[1]
            if bd <= 8:
                self.CH = 2.1 - 0.1 * bd
            else:self.CH = 1.3
            self.Fg = 1/2 * rho * self.Ug ** 2 * self.CH * z
        elif len(shape) == 3:
            # 箱型，参数为宽度高度和腹板倾角
            self.D = shape[1]
            bd = shape[0] / shape[1]
            if bd < 8:
                self.CH = 2.1 - 0.1 * bd
            else:self.CH = 1.3
            if shape[2] < 60:
                self.CH *= 1 - 0.005 * shape[2]
            else:self.CH *= 0.7
            self.Fg = 1/2 * rho * self.Ug ** 2 * self.CH * z
        else:
            # 桁架，参数为形状、实面积比、中心距、桁架高度、各杆件高度
            self.D = shape[4:]
            shape_l = np.array(shape[4:])
            if shape[0] == 1:
                # 形状为矩形或 H 形
                self.CH = np.interp(shape[1], l_CH['R'], l_CH['RorH'])
            elif shape[0] == 2:
                # 形状为圆形
                self.CH = []
                for i in shape_l:
                    du = i * self.Ud
                    if du > 6:
                        ch = np.interp(shape[1], l_CH['R'], l_CH['du>6'])
                        self.CH.append(ch)
                    else:
                        ch = np.interp(shape[1], l_CH['R'], l_CH['du<6'])
                        self.CH.append(ch)
                self.CH = np.array(self.CH)
            else:
                self.CH = 0
                print('请输入桁架形状，矩形或 H 型钢为 1， 圆形为 2')
            etas = [np.interp(shape[1], l_eta.iloc[0, 1:], l_eta.iloc[i, 1:]) for i in np.arange(1, 7)]
            eta = np.interp(shape[2]/shape[3], np.arange(1, 7), etas)
            self.CH *= eta
            self.Fg = 1/2 * rho * self.Ug ** 2 * self.CH * shape_l

    def ffr(self, cf=0.01, s=10):
        """计算主梁顺桥向荷载

        参数：
            cf：摩擦系数，查表 5.3.6 获取
            s：主梁周长(m)，桁架断面为梁体外轮廓周长
        """
        if self.l < 200:
            if self.shape > 3:
                self.Ffr = self.Fg * 0.5
            else: self.Ffr = self.Fg * 0.25
        else:
            self.Cf = cf
            self.s = s
            self.Ffr = 1/2 * self.rho * self.Ug ** 2 * cf * s


class PierWind():
    """桥墩风荷载

    参数：
        wind：设计风实例
        shape：[桥墩高度，断面形状，顺桥宽度，横桥宽度，*倒角半径，*架设情况]
            -单位：m
            -断面形状：直接取表 5.4.2-1 行数(矩形为前八行)
            -架设情况：0、1(未架设、已架设)
        rho：空气密度(kg/m3)
        ud：如果缺少设计风资料，用于估计风速值，默认地表为 A (m/s)
    """

    def __init__(self, wind, shape, rho=1.25, ud=25):
        self.wind = wind
        self.z = shape[0] * 0.65
        self.h = shape[0]
        self.rho = rho
        # 是否为直接提供风速值
        if wind == 0:
            self.Ud = ud
            self.GV = np.interp(shape[0], l_GV_h['H'], l_GV_h['A'])
        else:
            self.Ud = wind.kf * (self.z / 10) ** wind.alpha0 * wind.Us10
            self.GV = np.interp(shape[0], l_GV_h['H'], l_GV_h[wind.surface_class])
        # 计算静阵风
        self.Ug = self.GV * self.Ud

        if len(shape) == 6 and shape[5] == 1:
            hw = [40, 40]
        else:
            hw = [shape[0]/shape[2], shape[0]/shape[3]]
        
        w = np.array([shape[2], shape[3]])
        t = np.array([shape[3], shape[2]])
        self.CD = []
        # 获取阻力系数
        if shape[1] <= 8:
            for i in [0, 1]:
                cds = [np.interp(hw[i], l_CD.iloc[0, 1:], l_CD.iloc[j, 1:]) for j in np.arange(1, 10)]
                cd = np.interp(t[i]/w[i], l_CD.iloc[1:10, 0], cds)
                self.CD.append(cd)
        else:
            for i in [0, 1]:
                cd = np.interp(hw[i], l_CD.iloc[0, 1:], l_CD.iloc[shape[1]+1, 1:])
                self.CD.append(cd)
        
        self.CD = np.array(self.CD)
        # 考虑桥墩倒角
        if len(shape) >= 5:
            self.CD *= np.fmax(0.5, 1 - 1.5 * shape[4] / w)
        # 计算等效静阵风荷载
        self.Fg = 1 / 2 * self.rho * self.Ug ** 2 * self.CD * w
        

            













    
            
        