import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection
from sectionproperties.pre.pre import Material

# 定义一些常数
a152 = 140  # 15.2预应力束面积（mm）
p_con = 1395  # 张拉控制应力

as_p = 0.2  # 预应力束距上下缘距离
ap_p = 0.2  # 各层预应力束距离

as_p_h = 0.2  # 预应力钢束横向边缘距离
ap_p_h = 0.2  # 预应力束横向距离


# 定义钢束类


class Tendon:
    """定义某构件预应力钢束
    说明：截面以下缘作为基准面，即下缘坐标为 0，其余坐标为正值。

    参数：
    sec：截面对象
    top：单位(m),截面上缘坐标
    bot：单位(m),截面下缘坐标（默认为0）
    """

    def __init__(self, sec, top, bot=0):
        self.sec = sec  # 确定截面对象
        self.ac = sec.get_area()  # 面积
        self.center_y = sec.get_c()[1]  # 形心 y
        self.Ic = sec.get_ic()[0]  # 惯性矩（竖向）
        self.ic = sec.get_rc()[0]  # 回转半径（竖向）
        self.top = top  # 上缘坐标
        self.bot = bot  # 下缘坐标
        self.h = self.top - self.bot  # 截面高度
        self.yu = self.top - self.center_y  # 形心距上缘距离
        self.yb = self.center_y - self.bot  # 形心距下缘距离

        # 另外一些可能用到的性质
        self.side = 'b'  # 弯矩方向（默认为正弯矩）
        self.ep = np.linspace(-self.yb, self.yu, 51)  # 默认偏心距列表
        self.np = [1, 1, 1, 1]  # 默认预应力值列表
        self.h0 = self.h - 0.2  # 默认有效高度

    def plot_pro(self, ep_list, np_list, ep_test=None, np_test=None):
        """根据偏心距和预应力值进行绘图

        :param ep_list:偏心距列表
        :param np_list: 预应力值列表，应有四个列表
        :param ep_test: 偏心距测试值
        :param np_test: 预应力值测试值
        :return:
        """
        fig, ax = plt.subplots(figsize=(8, 8))
        plt.plot(ep_list, np_list[0], label='1')
        plt.plot(ep_list, np_list[1], label='2')
        plt.plot(ep_list, np_list[2], label='3')
        plt.plot(ep_list, np_list[3], label='4')
        # plt.plot(ep_list[1], np_list[1][0], '--', label='-1')
        # plt.plot(ep_list[1], np_list[1][1], '--', label='-2')
        # plt.plot(ep_list[1], np_list[1][2], '--', label='-3')
        # plt.plot(ep_list[1], np_list[1][3], '--', label='-4')

        if ep_test is None and np_test is None:
            pass
        else:
            plt.scatter(ep_test, np_test)

        np_max = np.fmax(np_list[0], np_list[1])
        np_max = np.fmax(np_max, np_list[2])
        plt.fill_between(
            ep_list,
            np_max,
            np_list[3],
            where=np_list[3] >= np_max,
            alpha=.25,
            interpolate=True)
        # ax.axvline(self.yb, ls=':', color='r')

        plt.ylim((0, max(np_list[3])))
        plt.xlabel('ep (m)')
        plt.ylabel('1/Np (1/kN)')
        plt.legend()
        plt.show()

    def tendon_cal(
            self,
            mg,
            me,
            ms,
            sig_t_1,
            sig_c_1,
            sig_t_2,
            sig_c_2,
            a=0.8,
            x=0.85):
        """获取最大弯矩对应截面所需预应力（正常使用极限状态）

        参数：
        Mg：单位(kN.m),恒载作用下弯矩列表 (正弯矩、负弯矩)
        Me：单位(kN.m),弹性组合下弯矩列表 (正弯矩、负弯矩)
        Ms：单位(kN.m),短期组合下弯矩列表 (正弯矩、负弯矩)
        sig_t_1：单位(MPa)，施工阶段混凝土容许拉应力
        sig_c_1：单位(MPa)，施工阶段容许压应力
        sig_t_2：单位(MPa)，使用阶段允许拉应力（抗裂）
        seg_c_2：单位(MPa)，使用阶段容许压应力
        a：预应力损失系数，默认 0.8
        x：预应力作用系数，默认 0.85（全预应力为 0.85，其余为 1）
        """
        sig_t_1 *= 1000
        sig_c_1 *= 1000
        sig_t_2 *= 1000
        sig_c_2 *= 1000
        # 确定正负弯矩

        if mg > 0:
            ep_list = -self.ep
            self.side = 'b'
            np_r1 = (ep_list * self.yu / (self.ic ** 2) - 1) / \
                    (self.ac * (mg * self.yu / self.Ic - sig_t_1))
            np_r2 = (1 + ep_list * self.yb / (self.ic ** 2)) / \
                    (self.ac * (mg * self.yb / self.Ic + sig_c_1))
            np_r3 = a * (1 - ep_list * self.yu / (self.ic ** 2)) / \
                (self.ac * (sig_c_2 - me / self.Ic * self.yu))
            np_r4 = x * a * (1 + ep_list * self.yb / (self.ic ** 2)) / \
                (self.ac * (ms / self.Ic * self.yb - sig_t_2))
        else:
            ep_list = self.ep
            self.side = 'u'
            np_r1 = (ep_list * self.yb / (self.ic ** 2) - 1) / \
                    (self.ac * (-mg * self.yb / self.Ic - sig_t_1))
            np_r2 = (1 + ep_list * self.yu / (self.ic ** 2)) / \
                    (self.ac * (-mg * self.yu / self.Ic + sig_c_1))
            np_r3 = a * (1 - ep_list * self.yb / (self.ic ** 2)) / \
                (self.ac * (sig_c_2 - (-me) / self.Ic * self.yb))
            np_r4 = x * a * (1 + ep_list * self.yu / (self.ic ** 2)) / \
                (self.ac * (-ms / self.Ic * self.yu - sig_t_2))

        self.np = [np_r1, np_r2, np_r3, np_r4]

        self.plot_pro(self.ep, self.np)

    def get_np_from_ep(self, ep_test):
        """
        通过偏心距得到预应力大小范围（内部工具函数）
        :param ep_test: 偏心距，单位 m，上部为正，下部为负
        :return: 预应力范围，单位 kN，[np_min, np_max]
        """

        ep_list = self.ep
        np_list = self.np
        np1_max_1 = np.interp(ep_test, ep_list, np_list[0])
        np1_max_2 = np.interp(ep_test, ep_list, np_list[1])
        np1_max_3 = np.interp(ep_test, ep_list, np_list[2])
        np_max = 1 / np.max([np1_max_1, np1_max_2, np1_max_3])
        np_min = 1 / np.interp(ep_test, ep_list, np_list[3])
        return [np_min, np_max]

    def get_np_from_pw(
            self,
            pro_width,
            s_p=as_p,
            p_p=ap_p,
            s_p_h=as_p_h,
            p_p_h=ap_p_h,
            pro_one=22):
        """通过预应力布置宽度等几何信息获得所需预应力

        参数：
        pw：单位(m),可布置预应力宽度
        a_p：单位(m),预应力距竖向边缘最小距离
        p_p：单位(m),预应力竖向最小距离
        ah_p：单位(m),预应力距横向边缘最小距离
        ph_p：单位(m),预应力横向最小距离
        pro_one：单束预应力股束，默认为 22 股
        :return: 预应力大小、预应力偏心距、预应力布置（从上到下）、预应力距边缘距离（从上到下）
        """
        pro_n = [0, ]
        pro_max_n = np.ceil((pro_width - 2 * s_p_h + 0.01) / p_p_h)

        bu = self.side
        if bu == 'b':
            print('对正弯矩(下缘)进行计算')
            pro_ep = [s_p - self.yb, ]
            edge = self.yb
            add_direction = 1
        elif bu == 'u':
            print('对正弯矩(上缘)进行计算')
            edge = self.yu
            pro_ep = [self.yu - s_p, ]
            add_direction = -1
        else:
            print('bu 变量只能为 b 或者 u，请重新输入')
            return 'wrong!'
        p_p = add_direction * p_p

        while True:
            ep_test = pro_ep[0]
            np_min = self.get_np_from_ep(ep_test)[0]
            np_max = self.get_np_from_ep(ep_test)[1]
            if np_max < np_min:
                print('起始偏心距过大!')
                pro_ep[0] += add_direction * 0.05
            elif pro_ep[0] * add_direction > 0:
                print('这个不行！')
                return 'not good!'
            else:
                print('起始偏心距合适...')
                break

        while True:
            if len(pro_n) > 10 or pro_one <= 0:
                print('无法完成配束')
                return

            np_test = sum(pro_n) * pro_one * a152 * p_con / 1000
            np_min = self.get_np_from_ep(ep_test)[0]
            np_max = self.get_np_from_ep(ep_test)[1]

            if np_min < np_test < np_max:
                break
            elif np_test < np_min:
                if pro_n[-1] == pro_max_n:
                    pro_n.append(1)
                    pro_ep.append(pro_ep[-1] + p_p)
                else:
                    pro_n[-1] += 1
                ep_test = np.average(pro_ep, weights=pro_n)
            elif np_test > np_max:
                pro_one -= 1

        self.plot_pro(self.ep, self.np, ep_test, 1 / np_test)

        pro_to_top = np.array(pro_ep) - self.yu
        print(f'预应力布置为：{pro_n}，合计{sum(pro_n)}束'
              f'每束{pro_one}根，{sum(pro_n) * pro_one}根')
        print(f'预应力总大小为：{np_test:.0f}kN')
        print(f'预应力整体偏心距为：{ep_test:.3f}m，具体距上缘距离为{pro_to_top}')

        pro_n = pro_n[::-add_direction]
        pro_to_top = pro_to_top[::-add_direction]

        return np_test, ep_test, pro_n, pro_to_top

    def get_np_from_mu(self, M, fc, b, ep, a1=1):
        """简单估算（不够准确）承载力所需钢束

        M：单位(kN.m)，基本组合弯矩
        fc：单位(Mpa)，混凝土抗压强度设计值
        b：单位(m)，截面宽度
        ep：单位(m)，预应力钢束偏心距
        a1：受压区系数，C50以下取 1.0，C50-C80 之间 1.0-0.94 内插

        :return：预应力大小，单位 kN
        """
        bu = self.side
        if bu == 'b':
            self.h0 = self.yu - ep
        elif bu == 'u':
            self.h0 = self.yb + ep
            M = -M

        np_out = b * 1000 * a1 * fc * \
                 (self.h0 * 1000 - ((self.h0 * 1000) ** 2 - 2 * M * 1000000 / (a1 * fc * b * 1000)) ** 0.5)
        print(f'需要预应力{np_out / 1000:.0f}kN')
        return np_out / 1000

    def get_ep_from_np(self, np_test, pro_n, s_p=as_p, p_p=ap_p):
        """通过预应力大小获得合适偏心距

        参数：
        np_test：单位(kN)，预应力大小
        pro_n：预应力束布置，从上到下
        """
        bu = self.side
        if bu == 'b':
            top_test = s_p - self.yb + p_p * (len(pro_n) - 1)
            edge = self.yb
            side = 1
        elif bu == 'u':
            top_test = self.yu - s_p
            edge = self.yu
            side = -1
        else:
            print('bu 变量只能为 b 或者 u，请重新输入')
            return 'wrong'

        top_good = []
        ep_good_list = []
        ep_good = []
        while True:
            ep_list = top_test - np.arange(len(pro_n)) * p_p
            ep_test = np.average(ep_list, weights=pro_n)
            np_min = self.get_np_from_ep(ep_test)[0]
            np_max = self.get_np_from_ep(ep_test)[1]

            if ep_test * side > 0:
                break
            elif np_min < np_test < np_max:
                top_good.append(top_test)
                ep_good_list.append(ep_list)
                ep_good.append(ep_test)

            top_test += side * 0.01

        if len(top_good) == 0:
            print('无合适偏心距满足要求！')
            return 'not good!'
        else:
            pro_to_top = np.array(top_good) - self.yu
            print(f'偏心距区间为：[{ep_good[0]:.2f}, {ep_good[-1]:.2f}]')
            print(f'上层钢束距上缘距离区间为[{pro_to_top[0]:.2f}, {pro_to_top[-1]:.2f}]')
            return ep_good[len(ep_good) // 2], ep_good_list[len(ep_good_list) // 2] - self.yu

# 测试
# beam_height = 2.8
# beam_width = 2.5
# g = sections.RectangularSection(d=beam_height, b=beam_width)
# g.plot_geometry()
# m = g.create_mesh(mesh_sizes=[0.01])
# s = CrossSection(g, m)
# s.plot_mesh()
# s.calculate_geometric_properties()
# s.calculate_warping_properties()
#
# test = Tendon(s, beam_height)
# mg = [31000, 46500]
# me = [36700, 54800]
# ms = [35600, 53100]
# m = [51154, 76219]
#
# # 1. 正常使用极限状态下预应力可选范围估算
# test.tendon_cal(mg, me, ms, 0, 22.41*1000, 0, 22.4*1000)
# # 2. 通过几何构造信息布置预应力（较大一边）
# pro_inf = test.get_np_from_pw(1.6, 'u')
