"""用于进行桩长计算

基于《公路桥涵地基与基础设计规范》(JTG D63-2007)进行计算"""

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']    # 绘图支持中文

import os

current_path = os.path.dirname(__file__)

# 一些常数


class Soil():
    """定义地基土的基本情况
    
    土体性质：置于 soil.csv 中，包括 name, height, depth, fa0, qik, frk, per
    rho：土层加权平均重度 (kN/m3)，默认为 18
    """

    def __init__(self, rho=18):
        """构造函数"""

        self.prop = pd.read_csv('soil.csv')
        self.rho = rho


class Pile_mc_zk():
    """定义摩擦钻孔桩基基本情况
    
    soil：土体对象
    d：桩基直径 (m)
    l：桩长 (m)
    t：桩端沉渣厚度 (m)，默认为 0.3
    h1：桩顶标高 (m)，默认为 0，注意与土体标高协调
    rho：桩基密度 (kN/m3)，默认为 26
    k2：地基土承载力深度修正系数，默认为 1.5
    """

    def __init__(self, soil, d, l, h1=0, rho=26, t=0, k2=1.5):
        """构造函数"""

        # 土体
        self.soil = soil
        # 桩基截面
        self.d = d
        self.u = d * np.pi
        self.ap = np.pi * (d / 2) ** 2
        # 桩长及密度
        self.l = l
        self.h1 = h1
        self.h2 = h1 - l
        self.rho = rho
        # 清底系数
        if t == 0:
            self.t = 0.3 * d
        self.m0 = np.interp(t/d, [0.3, 0.1], [0.7, 1])
        # 桩侧摩阻力
        self.ra1 = 0
        self.h0 = self.soil.prop['height'][0] + self.soil.prop['depth'][0]
        for i, j in self.soil.prop.iterrows():
            if j['height'] > self.h2:
                self.ra1 += 1/2 * self.u * j['qik'] * j['depth']
                self.permeable = j['per']
            else:
                self.ra1 += 1/2 * self.u * j['qik'] * (j['depth'] - (self.h2 - j['height']))
                self.fa0 = j['fa0']
                self.permeable = j['per']
                break
        # 修正系数
        if self.permeable:
            self.lam = np.interp(l/d, [20, 25], [0.7, 0.85])
        else:
            self.lam = np.interp(l/d, [20, 25], [0.65, 0.72])
        # 桩端承载力
        self.k2 = k2
        self.h_cal = min(self.h0 - self.h2, 40) - 3
        self.qr = self.m0 * self.lam * (self.fa0 + k2 * self.soil.rho * self.h_cal)
        self.ra2 = self.ap * self.qr
        # 承载力容许值
        self.ra0 = self.ra1 + self.ra2
        self.ra = self.ra0 - (self.ap * self.l * self.rho) + (self.h0 - self.h2) * self.ap * self.soil.rho


class Pile_mc_cz():
    """定义摩擦沉桩基基本情况
    
    soil：土体对象
    d：桩基直径 (m)
    l：桩基长度 (m)
    h1：桩顶标高 (m)，默认为 0
    rho：桩基密度 (kN/m3)，默认为 26
    alpha_i：振动沉桩对侧摩阻的影响系数，默认为 1
    alpha_r：振动沉桩对桩端承载力的影响系数，默认为 1
    """

    def __init__(self, soil, d, l, h1=0, rho=26, alpha_i=1, alpha_r=1):
        # 土体
        self.soil = soil
        # 桩基界面
        self.d = d
        self.u = d * np.pi
        self.ap = np.pi * (d / 2) ** 2
        # 桩长及密度
        self.l = l
        self.h1 = h1
        self.h2 = h1 - l
        self.rho = rho
        # 桩侧摩阻力
        self.alpha_i = alpha_i
        self.ra1 = 0
        self.h0 = self.soil.prop['height'][0] + self.soil.prop['depth'][0]
        for i, j in self.soil.prop.iterrows():
            if j['height'] > self.h2:
                self.ra1 += self.u * alpha_i * j['qik'] * j['depth']
            else:
                self.ra1 += self.u * alpha_i * j['qik'] * (j['depth'] - (self.h2 - j['height']))
                self.qrk = j['fa0']
        # 桩端承载力
        self.alpha_r = alpha_r
        self.ra2 = self.alpha_r * self.ap * self.qrk
        # 承载力容许值
        self.ra1 *= 1 / 2
        self.ra2 *= 1 / 2
        self.ra0 = self.ra1 + self.ra2
        self.ra = self.ra0 - (self.ap * self.l * self.rho) + (self.h0 - self.h2) * self.ap * self.soil.rho


class Pile_dc():
    """定义端承桩基基本情况
    
    soil：土体对象
    d：桩基直径 (m)
    l：桩基长度 (m)
    h1：桩顶标高 (m)，默认为 0
    rho：桩基密度 (kN/m3)，默认为 26
    type：桩基类型，1-钻孔桩、2-沉桩，默认为钻孔桩
    com：岩石层情况，1-完整、2-较破碎、3-破碎，默认为 1；
    """

    def __init__(self, soil, d, l, h1=0, rho=26, type=1, com=1):
        # 土体
        self.soil = soil
        # 桩基界面
        self.d = d
        self.u = d * np.pi
        self.ap = np.pi * (d / 2) ** 2
        # 桩长及密度
        self.l = l
        self.h1 = h1
        self.h2 = h1 - l
        self.rho = rho
        # c1/c2 系数取值
        self.completion = com
        if com == 1:
            self.cc1 = 0.6
            self.cc2 = 0.05
        elif com == 2:
            self.cc1 = 0.5
            self.cc2 = 0.04
        elif com == 3:
            self.cc1 = 0.4
            self.cc2 = 0.03
        else:
            print('输入错误！')
        # 根据桩基类型折减 c1/c2
        self.type = type
        if type == 1:
            self.cc1 *= 0.8
            self.cc2 *= 0.8
        self.c1 = []
        self.c2 = []
        
        # 各分项力
        self.ra1 = 0    # 土层摩阻力
        self.ra2 = 0    # 岩层摩阻力
        self.ra3 = 0    # 桩端承载力
        self.h0 = self.soil.prop['height'][0] + self.soil.prop['depth'][0]

        for i, j in self.soil.prop.iterrows():
            if j['height'] > self.h2:
                self.c1.append(self.cc1)
                self.c2.append(self.cc2)
                if j['frk'] > 0:
                    if j['depth'] <= 0.5: self.c2[-1] = 0
                    if '中风化' in j['name']: self.c2[-1] *= 0.75

                    self.ra2 += self.c2[-1] * self.u * j['frk'] * j['depth']
                else:
                    self.ra1 += 1/2 * self.u * j['qik'] * j['depth']
            else:
                self.c1.append(self.cc1)
                self.c2.append(self.cc2)
                if j['frk'] >= 2000:
                    if j['frk'] < 15000: self.xi = 0.8
                    elif j['frk'] < 30000: self.xi = 0.5
                    else: self.xi = 0.2
                    
                    if j['depth'] - (self.h2 - j['height']) < 0.5:
                        self.c1[-1] *= 0.75
                        self.c2[-1] = 0
                    if '中风化' in j['name']:
                        self.c1[-1] *= 0.75
                        self.c2[-1] *= 0.75
                    
                    self.ra2 +=  self.c2[-1] * self.u * j['frk'] * (j['depth'] - (self.h2 - j['height']))
                    self.ra3 += self.c1[-1] * self.ap * j['frk']
                else:
                    self.xi = 1
                    self.ra1 += 1/2 * self.u * j['qik'] * (j['depth'] - (self.h2 - j['height']))
                break
        # 端部承载力
        self.ra1 *= self.xi
        self.ra0 = self.ra1 + self.ra2 + self.ra3
        self.ra = self.ra0 - (self.ap * self.l * self.rho) + (self.h0 - self.h2) * self.ap * self.soil.rho


def get_l(soil, d, F, factor=1.25, h1=0, rho=26, t=0, k2=1.5, type=1, com=1):
    """用于根据安全系数确定桩长

    soil：土体对象
    d：桩基直径 (m)
    F：桩顶力 (kN)
    factor：安全系数
    """

    l = [0,]
    ra_1 = [0,]
    ra_2 = [0,]
    while ra_1[-1] <= F * factor and ra_2[-1] <= F * factor:
        l.append(l[-1] + 0.5)
        pile_1 = Pile_mc_zk(soil, d, l[-1], h1, rho, t, k2)
        pile_2 = Pile_dc(soil, d, l[-1], h1, rho, type, com)
        ra_1.append(pile_1.ra)
        ra_2.append(pile_2.ra)
        if pile_1.h2 - 0.5 < pile_1.soil.prop.height.min():
            print('土层可能太浅')
            break
    return np.array(l)


def pile_l(llist, soil, d, F, h1=0, rho=26, t=0, k2=1.5, type=1, com=1):
    """ 用于绘制桩基承载力随桩长的变化关系

    llist：绘制桩长列表
    soil：土体对象
    d：桩基直径
    F：桩顶力
    """

    pile_ra_1 = []
    pile_ra_2 = []
    for i in llist:
        pile_1 = Pile_mc_zk(soil, d, i, h1, rho, t, k2)
        pile_2 = Pile_dc(soil, d, i, h1, rho, type, com)
        pile_ra_1.append(pile_1.ra)
        pile_ra_2.append(pile_2.ra)
    pile_ra_1 = np.array(pile_ra_1)
    pile_ra_2 = np.array(pile_ra_2)
    
    plt.figure(figsize=(10, 3))

    ax = plt.subplot(1, 2, 1)
    ax.set_yticks([F, ], minor = True)
    ax.yaxis.grid(True, which ='minor', ls='--', c='r')
    # plt.plot(llist, F * np.ones_like(llist), 'r--',linewidth=0.7, label='桩顶力')
    plt.plot(llist, pile_ra_1, label='摩擦桩')
    plt.plot(llist, pile_ra_2, label='端承桩')
    plt.legend(loc='lower right')
    plt.xlabel('桩长 (m)')
    plt.ylabel('单桩承载力 (kN)')
    plt.text(llist[-1], pile_ra_1[-1], f'({llist[-1]:.1f}, {pile_ra_1[-1]:.0f})', ha='right', va='center')
    plt.text(llist[-1], pile_ra_2[-1], f'({llist[-1]:.1f}, {pile_ra_2[-1]:.0f})', ha='right', va='center')
    plt.text(llist[0], F*1.05, f'桩顶力：{F:.0f}', ha='left')

    ax2 = plt.subplot(1, 2, 2)
    ax2.set_yticks([1.001, ], minor = True)
    ax2.yaxis.grid(True, which ='minor', ls='--', c='r')
    plt.plot(llist, pile_ra_1/F, label='摩擦桩')
    plt.plot(llist, pile_ra_2/F, label='端承桩')
    plt.legend(loc='lower right')
    plt.xlabel('桩长 (m)')
    plt.ylabel('安全系数')
    plt.text(llist[-1], pile_ra_1[-1]/F, f'({llist[-1]:.1f}, {pile_ra_1[-1]/F:.2f})', ha='right', va='center')
    plt.text(llist[-1], pile_ra_2[-1]/F, f'({llist[-1]:.1f}, {pile_ra_2[-1]/F:.2f})', ha='right', va='center')
    return pile_ra_1, pile_ra_2

