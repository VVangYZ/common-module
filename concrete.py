import sympy
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def get_as(rebar_d, num):
    return np.pi * (rebar_d / 2) ** 2 * num


class RectangleCompress:
    def __init__(self, b, h, ad1, ad2, anum1, anum2, fcd, fsd1, fsd2, a1, a2):
        """
        矩形截面偏心受压构件计算
        :param b: 矩形宽度，m
        :param h: 矩形高度，m
        :param ad1: 受拉钢筋直径，mm
        :param ad2: 受压钢筋直径，mm
        :param anum1: 受拉钢筋根数
        :param anum2: 受压钢筋根数
        :param fcd: 混凝土抗压强度设计值，Mpa
        :param fsd1: 钢筋抗拉强度设计值，Mpa
        :param fsd2: 钢筋抗压强度设计值，MPa
        :param a1: 受拉钢筋中心到边缘距离，m
        :param a2: 受压钢筋中心到边缘距离，m
        """
        self.b = b
        self.h = h
        self.ad1 = ad1
        self.ad2 = ad2
        self.anum1 = anum1
        self.anum2 = anum2
        self.as1 = get_as(ad1, anum1)
        self.as2 = get_as(ad2, anum2)
        self.a1 = a1
        self.a2 = a2
        self.fcd = fcd
        self.fsd1 = fsd1
        self.fsd2 = fsd2

        self.xi_b = 0.49  # 相对界限受压区高度
        self.h0 = h - a1  # 截面有效高度

    # noinspection PyAttributeOutsideInit
    def capacity(self, nd, md, r=1.1):
        self.e0 = md / nd  # 截面偏心距
        self.e = self.e0 + self.h / 2 - self.a1  # 截面相对受拉钢筋偏心距
        self.m_load = r * nd * self.e

        def nud_cal(x):
            nud = self.fcd * self.b * x * 1e3 + self.fsd2 * self.as2 / 1e3 - self.fsd1 * self.as1 / 1e3
            return nud

        def mud_cal(x):
            mud = self.fcd * self.b * x * (self.h0 - x / 2) * 1e3
            mud += self.fsd2 * self.as2 * (self.h0 - self.a2) / 1e3
            return mud

        def e0_cal(x):
            return mud_cal(x) - self.e * nud_cal(x)

        self.x = optimize.root(e0_cal, self.h0 * 0.5).x[0]

        # x = sympy.Symbol('x')
        # self.x = sympy.solve(r * nd * 1e3 - self.fcd * self.b * x * 1e6 -
        #                      self.fsd2 * self.as2 + self.fsd1 * self.as1, x)[0]

        self.xi = self.x / self.h0  # 截面相对受压区高度
        self.type = 'Large' if self.xi <= self.xi_b else 'Small'
        self.n_resistance = nud_cal(self.x)
        self.m_resistance = mud_cal(self.x)
        # self.m_resistance = self.fcd * self.b * self.x * (self.h0 - self.x / 2) * 1e3 + \
        #                     self.fsd2 * self.as2 * (self.h0 - self.a2) / 1e3

        return self.m_load, self.m_resistance

    # noinspection PyAttributeOutsideInit
    def crack_width(self, ns, ms, es=2e5, rebar=36):
        c1 = 1
        c2 = 1.5 * 0.9
        c3 = 0.9
        self.e0s = ms / ns  # 频遇组合下偏心距

        if self.e0s <= 0.55 * self.h:
            self.wcr = 0
        else:
            self.es = self.e0s + self.h / 2 - self.a1  # 频遇组合下
            z = (0.87 - 0.12 * (self.h0 / self.es) ** 2) * self.h0
            self.sigma_ss = ns * (self.es - z) / (self.as1 * z) * 1000
            pte = min(self.as1 / (2 * self.a1 * self.b) / 1e6, 0.1)
            self.wcr = c1 * c2 * c3 * self.sigma_ss / es * ((self.a1 - rebar / 2000) + rebar) / (0.3 + 1.4 * pte)
        return self.wcr


class CircularCompress:
    def __init__(self, d, anum, ad, fcd, fsd1, fsd2, c=50, hoop=12):
        """
        定义圆形截面偏压构件
        :param d: 桩基直径（m）
        :param anum: 纵向钢筋数量
        :param ad: 纵向钢筋直径（mm）
        :param fcd: 混凝土抗压强度设计值（MPa）
        :param fsd1: 钢筋抗拉强度设计值（MPa）
        :param fsd2: 钢筋抗压强度设计值（MPa）
        :param c: 保护层厚度 （mm）
        :param hoop：箍筋直径 （mm）
        """
        self.d = d
        self.r = d / 2
        self.a = np.pi * (d / 2) ** 2
        self.anum = anum
        self.ad = ad
        self.As = anum * np.pi * (self.ad / 2) ** 2
        self.fcd = fcd
        self.fsd1 = fsd1
        self.fsd2 = fsd2
        self.Es = 2e5

        self.c = c
        self.hoop = hoop
        self.cs = c + hoop
        self.rs = self.r - (self.cs + self.ad / 2) / 1000

        self.e0 = 0
        self.cap_axes = 0
        self.nud = 0
        self.mud = 0
        self.alpha = 0

        self.wcr = 0

        self.vr = 0

    def capacity(self, nd, md, r=1.1):
        """
        计算截面承载力
        :param nd: 基本组合下轴力（kN）
        :param md: 基本组合下弯矩（kN.m)
        :param r: 结构重要性系数，默认 1.1
        :return:
        """
        self.cap_axes = 0.9 * (self.fcd * self.a * 1e6 + self.As * self.fsd2) / 1000
        if self.cap_axes < r * nd:
            print('轴心抗压承载力无法通过！')
            raise Exception

        self.e0 = md / nd

        # alpha = 0
        # alpha_t = 0

        def nud_cal(alpha):
            alpha_t = 0 if alpha > 0.625 else (1.25 - 2 * alpha)
            nud = alpha * self.fcd * self.a * (1 - np.sin(2 * np.pi * alpha) / (2 * np.pi * alpha)) * 1e3
            nud += (alpha - alpha_t) * self.fsd1 * self.As / 1e3
            return nud
        
        def mud_cal(alpha):
            alpha_t = 0 if alpha > 0.625 else (1.25 - 2 * alpha)
            mud = 2 / 3 * self.fcd * self.a * self.d / 2 * (np.sin(alpha * np.pi) ** 3) / np.pi * 1e3
            mud += self.fsd1 * self.As * self.rs / 2 * (np.sin(alpha * np.pi) + np.sin(np.pi * alpha_t)) / np.pi / 1e3
            return mud
        
        def e0_cal(alpha):
            return mud_cal(alpha) - self.e0 * nud_cal(alpha)
        
        alpha = optimize.root(e0_cal, 0.5).x[0]
        self.alpha = alpha
        self.nud = nud_cal(alpha)
        self.mud = mud_cal(alpha)
        if self.nud < r * nd:
            print('计算不通过！')
            print(alpha)

    def crack_width(self, ns, ms):
        """

        :param ns: 短期组合轴力（kN)
        :param ms: 短期组合弯矩 (kN.m)
        :return: 裂缝计算宽度 (mm)
        """
        c1 = 1
        c2 = 1.5 * 0.9
        c3 = 0.9

        c = self.cs       # 保护层厚度加箍筋直径
        es = self.r - self.rs        # 纵向钢筋到边缘距离
        rs = self.r - es     # 钢筋布置半径
        e0 = ms / ns        # 截面偏心距

        if e0 <= 0.55 * self.r:     # 是否需要验算裂缝
            self.wcr = 0
            return

        rho = self.As / 1e6 / (np.pi * self.r ** 2)     # 配筋率
        beta = (0.4 + 2.5 * rho) * (1 + 0.353 * (e0 / self.r) ** -2)
        r1 = self.r - 2 * es
        rho_te = beta * self.As / 1e6 / (np.pi * (self.r ** 2 - r1 ** 2))     # 有效配筋率

        n_s = 1 + 1 / (4000 * e0 / (2 * self.r - es / 1e3)) * (20 / self.d) ** 2
        sigma_ss = 0.6 * (e0 / self.r - 0.1) ** 3 / ((0.45 + 0.26 * rs / self.r) * (n_s * e0 / self.r + 0.2) ** 2)
        sigma_ss *= ns * 1e3 / self.As

        c = min(c, 50)
        self.wcr = c1 * c2 * c3 * sigma_ss / self.Es * (c + self.ad / 2) / (0.36 + 1.7 * rho_te)

        return self.wcr

    def shear_capacity(self, fc1, sk, vc0, fyh=0):
        """
        计算抗剪承载力
        :param fc1: 混凝土抗压强度标准值（MPa）
        :param sk: 箍筋间距（mm）
        :param vc0：剪力设计值（kN）
        :param fyh：箍筋抗拉强度设计值（MPa，默认为主筋强度）
        :return:
        """
        fyh = fyh if fyh != 0 else self.fsd1
        de = self.d - (self.c + 0.5 * self.hoop) / 1e3      # 核心混凝土直径（m）
        ae = (de / 2) ** 2 * np.pi * 1e4                    # 核心混凝土面积（cm2）
        vc = 0.0023 * fc1 ** 0.5 * ae                       # 混凝土抗力（kN）
        ak = (self.hoop / 2) ** 2 * np.pi * 2 / 1e2         # 同一截面上箍筋面积（cm2）
        vs = min(0.1 * ak * (0.9 * self.d * 1e2) * fyh / (sk / 10), 0.067 * fc1 ** 0.5 * ae)
        self.vr = vc + vs
        safe = abs(self.vr / vc0)

        return safe


class RectangleBend:
    def __init__(self, b, h, fcd, ad1, anum1, fsd1, ad2, anum2, fsd2, pd=15.4, pnum=0, fpd=1260, a=50, hoop=16):
        """
        :param b: 矩形截面宽度，m
        :param h: 矩形截面高度，m
        :param fcd: 混凝土抗压强度设计值，MPa

        :param ad1: 受拉钢筋直径，mm
        :param anum1: 受拉钢筋根数
        :param fsd1: 钢筋抗拉强度设计值，MPa
        :param ad2: 受压钢筋直径，mm
        :param anum2: 受压钢筋根数
        :param fsd2: 钢筋抗压强度设计值，MPa
        :param pd: 预应力钢束直径，mm
        :param pnum: 预应力钢束根数
        :param fpd: 预应力钢束抗拉强度设计值，MPa
        :param a: 结构保护层厚度，mm
        :param hoop: 箍筋直径，mm
        """
        self.b = b
        self.h = h
        self.ad1 = ad1
        self.ad2 = ad2
        self.pd = pd
        self.anum1 = anum1
        self.anum2 = anum2
        self.pnum = pnum
        self.as1 = get_as(ad1, anum1)
        self.as2 = get_as(ad2, anum2)
        self.ap = get_as(pd, pnum)
        self.a = a
        self.fcd = fcd
        self.fsd1 = fsd1
        self.fsd2 = fsd2
        self.fpd = fpd
        self.hoop = hoop

        self.xi_b = 0.49
        self.aa = a + hoop + ad1 / 2        # 钢筋合力点到边缘
        self.h0 = h - self.aa / 1e3

        self.x = 0
        self.m_resistance = 0

        self.vcs = 0
        self.vpd = 0
        self.vr = 0

    def capacity(self, md=1):
        """

        :param md: 弯矩设计值，kN.m
        :return: 安全系数
        """
        self.x = (self.fsd1 * self.as1 + self.fpd * self.ap - self.fsd2 * self.as2) / 1e6 / (self.fcd * self.b)
        if self.x > self.xi_b * self.h0:
            print('超筋破坏')
            raise Exception
        elif self.x < 2 * self.aa / 1e3:
            print('受压钢筋未屈服')
            self.m_resistance = self.fpd * self.ap * (self.h0 - 0.25) / 1e3
            self.m_resistance += self.fsd1 * self.as1 * (self.h0 - self.aa / 1e3) / 1e3
        else:
            self.m_resistance = self.fcd * self.b * self.x * (self.h0 - self.x) * 1e3
            self.m_resistance += self.fsd2 * self.as2 * (self.h0 - self.aa / 1e3) / 1e3
        return self.m_resistance / md

    def shear_capacity(self, vd=1, space_sv=100, fcuk=50, fsv=330, num_sv=1, p_theta=0):
        psv = get_as(self.hoop, 2) * num_sv / (space_sv * self.b * 1e3)
        a_2 = 1.25 if self.pnum != 0 else 1
        p = (self.as1 + self.ap) / (self.b * self.h0) / 1e6 * 100
        p = min(p, 2.5)
        self.vcs = 0.45 * 1e-3 * a_2 * self.b * self.h0 * 1e6 * ((2 + 0.6 * p) * fcuk ** 0.5 * (psv * fsv)) ** 0.5
        self.vpd = 0.75 * 1e-3 * self.fpd * self.ap * np.sin(p_theta * np.pi / 180)
        self.vr = self.vcs + self.vpd
        return self.vr / vd


if __name__ == '__main__':
    # a = CircularCompress(1.6, 29, 32, 13.8, 415, 400)
    # a.capacity(3800, 5800)
    # a.nud / (3800 * 1.1)
    # a.crack_width(6000, 5000)
    # a.shear_capacity(30, 150, 1000)
    # b = RectangleCompress(1.8, 1.8, 32, 32, 17, 17, 22.4, 415, 400, 0.094, 0.094)
    # b.capacity(3434, 7867, 1)
    # print(b.m_load, b.m_resistance)
    c = RectangleBend(
        b=2.1,
        h=2.2,
        fcd=22.4,
        ad1=32,
        anum1=14,
        fsd1=415,
        ad2=32,
        anum2=16,
        fsd2=400,
        a=40,
        pnum=10*15
    )
    c.capacity(10000)
    c.shear_capacity(num_sv=1.5 + 3 / 4, p_theta=8.7)
    c.m_resistance
    c.vr




