import numpy as np
import pandas as  pd
import matplotlib.pyplot as plt

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection
from sectionproperties.pre.pre import Material

# 定义一些常数
a152 = 140      # 15.2预应力束面积（mm）
pcon = 1395       # 张拉控制应力考虑 20% 预应力损失

asp = 0.2       # 预应力束距上下缘距离
aap = 0.2       # 各层预应力束距离

asph = 0.2      # 预应力钢束横向边缘距离
aaph = 0.2      # 预应力束横向距离

# 定义钢束类
class Tendon(): 
    """定义某构件预应力钢束

    参数：
    sec：截面对象
    top：单位(m),截面上缘坐标
    bot：单位(m),截面下缘坐标（默认为0）
    """

    def __init__(self, sec, top, bot=0):
        self.sec = sec
        self.ac = sec.get_area()        # 面积
        self.center_y = sec.get_c()[1]      # 形心 y
        self.Ic = sec.get_ic()[0]       # 惯性矩
        self.ic = sec.get_rc()[0]       # 回转半径
        self.top = top      # 上缘坐标
        self.bot = bot      # 下缘坐标
        self.yu = self.top - self.center_y      # 形心距上缘距离
        self.yb = self.center_y - self.bot      # 形心距下缘距离


    def plot_pro(self, ep_list, np_list, ep_test=0, np_test=0):
        fig = plt.subplots(figsize=(8, 8))
        plt.plot(ep_list[0], np_list[0][0], label='1')
        plt.plot(ep_list[0], np_list[0][1], label='2')
        plt.plot(ep_list[0], np_list[0][2], label='3')
        plt.plot(ep_list[0], np_list[0][3], label='4')
        plt.plot(ep_list[1], np_list[1][0], '--', label='-1')
        plt.plot(ep_list[1], np_list[1][1], '--', label='-2')
        plt.plot(ep_list[1], np_list[1][2], '--', label='-3')
        plt.plot(ep_list[1], np_list[1][3], '--', label='-4')

        if ep_test == 0 and np_test == 0:
            pass
        else:
            plt.scatter(ep_test, np_test)

        np_rmax_b = np.fmax(np_list[0][0], np_list[0][1])
        np_rmax_b = np.fmax(np_rmax_b, np_list[0][2])
        plt.fill_between(ep_list[0], np_rmax_b, np_list[0][3], where=np_list[0][3] >= np_rmax_b, alpha=.25, interpolate=True)
        # ax.axvline(self.yb, ls=':', color='r')
        np_rmax_u = np.fmax(np_list[1][0], np_list[1][1])
        np_rmax_u = np.fmax(np_rmax_u, np_list[1][2])
        plt.fill_between(ep_list[1], np_rmax_u, np_list[1][3], where=np_list[1][3] >= np_rmax_u, alpha=.25, interpolate=True)

        plt.xlabel('ep')
        plt.ylabel('1/Np')
        plt.legend()
        plt.show()


    def tendon_cal(self, Mg, Me, Ms, sigt1, sigc1, sigt2, sigc2, a=0.8, x=0.85):
        """获取最大弯矩对应截面所需预应力（正常使用极限状态）

        参数：
        Mg：单位(kN.m),恒载作用下弯矩列表 (正弯矩、负弯矩)
        Me：单位(kN.m),弹性组合下弯矩列表 (正弯矩、负弯矩)
        Ms：单位(kN.m),短期组合下弯矩列表 (正弯矩、负弯矩)
        sigt1：单位(kPa)，施工阶段混凝土容许拉应力
        sigc1：单位(kPa)，施工阶段容许压应力
        sigt2：单位(kPa)，使用阶段允许拉应力（抗裂）
        segc2：单位(kPa)，使用阶段容许压应力
        a：预应力损失系数，默认 0.8
        x：预应力作用系数，默认 0.85（全预应力为 0.85，其余为 1）
        """
        # 正弯矩（跨中）计算
        ep_test_b = np.linspace(0, self.yb, 31)
        np_r1 = (ep_test_b * self.yu / (self.ic ** 2) - 1) / (self.ac * (Mg[0] * self.yu/self.Ic - sigt1))
        np_r2 = (1 + ep_test_b * self.yb / (self.ic ** 2)) / (self.ac * (Mg[0] * self.yb/self.Ic + sigc1))
        np_r3 = a * (1 - ep_test_b * self.yu / (self.ic ** 2)) / (self.ac * (sigc2 - Me[0] / self.Ic * self.yu))
        np_r4 = x * a * (1 + ep_test_b * self.yb / (self.ic ** 2)) / (self.ac * Ms[0] / self.Ic * self.yb)

        # 负弯矩（支点）计算
        ep_test_u = np.linspace(0, self.yu, 31)
        np_r5 = (ep_test_u * self.yb / (self.ic ** 2) - 1) / (self.ac * (Mg[1] * self.yb / self.Ic - sigt1))
        np_r6 = (1 + ep_test_u * self.yu / (self.ic ** 2)) / (self.ac * (Mg[1] * self.yu / self.Ic + sigc1))
        np_r7 = a * (1 - ep_test_u * self.yb / (self.ic ** 2)) / (self.ac * (sigc2 - Me[1] / self.Ic * self.yb))
        np_r8 = x * a * (1 + ep_test_u * self.yu / (self.ic ** 2)) / (self.ac * Ms[1] / self.Ic * self.yu)

        self.ep_b = ep_test_b
        self.ep_u = ep_test_u
        self.np_b = [np_r1, np_r2, np_r3, np_r4]
        self.np_u = [np_r5, np_r6, np_r7, np_r8]
        self.big_side = 'b' if Ms[0] > Ms[1] else 'u'
        self.small_side = 'b' if Ms[0] < Ms[1] else 'u'

        self.plot_pro([self.ep_b, self.ep_u], [self.np_b, self.np_u])
        # self.np_b = [np_r4, np_rmax_b]
        # self.np_u = [np_r8, np_rmax_u]


    def get_np_from_ep(self, ep_test):
        """
        通过偏心距得到预应力范围
        :param ep_test: 偏心距，单位 m，上部为正，下部为负
        :return: 预应力范围，单位 kN，[np_min, np_max]
        """

        if ep_test > 0:
            ep_list = self.ep_u
            np_list = self.np_u
        else:
            ep_list = self.ep_b
            np_list = self.np_b
        np1_max_1 = np.interp(abs(ep_test), ep_list, np_list[0])
        np1_max_2 = np.interp(abs(ep_test), ep_list, np_list[1])
        np1_max_3 = np.interp(abs(ep_test), ep_list, np_list[2])
        np_max = 1 / np.max([np1_max_1, np1_max_2, np1_max_3])
        np_min = 1 / np.interp(abs(ep_test), ep_list, np_list[3])
        return [np_min, np_max]


    def get_np_from_pw(self, pw, bu=0, a_p=asp, p_p=aap, ah_p=asph, ph_p=aaph, pro_pern=22):
        """通过预应力布置宽度等几何信息获得所需预应力

        参数：
        pw：单位(m),可布置预应力宽度
        bu：顶部或者底部，bu='b'则为底部，bu='u'则为顶部，默认为 big_side
        a_p：单位(m),预应力距竖向边缘最小距离
        p_p：单位(m),预应力竖向最小距离
        ah_p：单位(m),预应力距横向边缘最小距离
        ph_p：单位(m),预应力横向最小距离
        pro_pern：单束预应力股束，默认为 22 股
        :return: 预应力大小、预应力偏心距、预应力布置（从上到下）、预应力距边缘距离（从上到下）
        """
        pro_n = [0, ]
        pro_max_n = np.ceil((pw - 2*ah_p) / ph_p)

        bu = self.big_side if bu == 0 else bu
        if bu == 'b':
            print('对正弯矩(下缘)进行计算')
            pro_ep = [a_p - self.yb,]
            edge = self.yb
            side = 1
        elif bu == 'u':
            print('对正弯矩(上缘)进行计算')
            edge = self.yu
            pro_ep = [self.yu - a_p,]
            side = -1
        else:
            print('bu 变量只能为 b 或者 u，请重新输入')
            return 'wrong!'
        p_p = side * p_p

        while True:
            ep_test = pro_ep[0]
            np_min = self.get_np_from_ep(ep_test)[0]
            np_max = self.get_np_from_ep(ep_test)[1]
            if np_max < np_min:
                print('起始偏心距过大!')
                pro_ep[0] += side * 0.05
            elif pro_ep[0] * side > 0:
                print('这个不行！')
                return 'not good!'
            else:
                print('起始偏心距合适...')
                break

        while True:
            if len(pro_n) > 10 or pro_pern <= 0:
                print('无法完成配束')
                return

            np_test = sum(pro_n) * pro_pern * a152 * pcon / 1000
            np_min = self.get_np_from_ep(ep_test)[0]
            np_max = self.get_np_from_ep(ep_test)[1]

            if np_test > np_min and np_test < np_max:
                break
            elif np_test < np_min:
                if pro_n[-1] == pro_max_n:
                    pro_n.append(1)
                    pro_ep.append(pro_ep[-1] + p_p)
                else:
                    pro_n[-1] += 1
                ep_test = np.average(pro_ep, weights=pro_n)
            elif np_test > np_max:
                pro_pern -= 1

        self.plot_pro([self.ep_b, self.ep_u], [self.np_b, self.np_u], abs(ep_test), 1 / np_test)

        pro_a_edge = np.array(pro_ep) - self.yu
        print(f'预应力布置为：{pro_n}，合计{sum(pro_n)}束'
              f'每束{pro_pern}根，{sum(pro_n) * pro_pern}根')
        print(f'预应力总大小为：{np_test:.0f}kN')
        print(f'预应力整体偏心距为：{ep_test:.3f}m，具体距上缘距离为{pro_a_edge}')

        pro_n = pro_n[::-side]
        pro_a_edge = pro_a_edge[::-side]

        return np_test, ep_test, pro_n, pro_a_edge


    def get_np_from_mu(self, M, fc, b, ep, bu=0, a1=1):
        """简单估算（不够准确）承载力所需钢束

        M：单位(kN.m)，基本组合弯矩
        fc：单位(Mpa)，混凝土抗压强度设计值
        b：单位(m)，截面宽度
        ep：单位(m)，预应力钢束偏心距
        bu：顶部或者底部，bu='b'则为底部，bu='u'则为顶部，默认为 big_side
        a1：受压区系数，C50以下取 1.0，C50-C80 之间 1.0-0.94 内插

        :return：预应力大小，单位 kN
        """
        bu = self.big_side if bu == 0 else bu
        if bu == 'b':
            self.h0 = self.yu - ep
        elif bu == 'u':
            self.h0 = self.yb + ep
        self.b = b
        np_out = b*1000 * a1 * fc * (self.h0*1000 - ((self.h0*1000)**2 - 2*M*1000000/(a1*fc*b*1000)) ** 0.5)
        print(f'需要预应力{np_out / 1000:.0f}kN')
        return np_out / 1000
   

    def get_ep_from_np(self, np_test, pro_n, bu=0, a_p=asp, p_p=aap):
        """通过预应力束数（15.2mm）获得合适偏心距

        参数：
        np_test：单位(kN)，预应力大小
        pro_n：预应力束布置，从上到下
        bu：顶部或者底部，bu='b'则为底部，bu='u'则为顶部，默认为 small_side
        """
        bu = self.small_side if bu == 0 else bu
        if bu == 'b':
            top_test =  a_p - self.yb + p_p * (len(pro_n) - 1)
            edge = self.yb
            side = 1
        elif bu == 'u':
            top_test = self.yu - a_p
            edge = self.yu
            side = -1
        else:
            print('bu 变量只能为 b 或者 u，请重新输入')

        top_good = []
        ep_goodlist = []
        ep_good = []
        while True:
            ep_list = top_test - np.arange(len(pro_n)) * p_p
            ep_test = np.average(ep_list, weights=pro_n)
            np_min = self.get_np_from_ep(ep_test)[0]
            np_max = self.get_np_from_ep(ep_test)[1]

            if ep_test * side > 0:
                break
            elif np_test > np_min and np_test < np_max:
                top_good.append(top_test)
                ep_goodlist.append(ep_list)
                ep_good.append(ep_test)

            top_test += side * 0.01

        if len(top_good) == 0:
            print('无合适偏心距满足要求！')
            return 'not good!'
        else:
            pro_to_top = np.array(top_good) - self.yu
            print(f'偏心距区间为：[{ep_good[0]:.2f}, {ep_good[-1]:.2f}]')
            print(f'上层钢束距上缘距离区间为[{pro_to_top[0]:.2f}, {pro_to_top[-1]:.2f}]')
            return  ep_goodlist[len(ep_goodlist) // 2] - self.yu


    def get_Np_from_ap(self, ap_b_test, ap_u_test, apro=a152, scon=pcon):
        """通过距离上下缘距离获得合适预应力大小

        参数：
        ap_b_test：单位(m),测试预应力作用点距下缘距离
        ap_u_test：单位(m),测试预应力作用点距上缘距离
        apro：单位（mm2），默认为 15.2 钢绞线截面积 140mm2
        scon：预应力张拉力，默认为 1395MPa
        """
        pro_N_1 = apro * scon / 1000
        ep_b_test = self.yb - ap_b_test
        ep_u_test = self.yu - ap_u_test
        np_b_min = 1 / np.interp(ep_b_test, self.ep_b, self.np_b[0])
        np_b_max = 1 / np.interp(ep_b_test, self.ep_b, self.np_b[1])
        
        np_u_min = 1 / np.interp(ep_u_test, self.ep_u, self.np_u[0])
        np_u_max = 1 / np.interp(ep_u_test, self.ep_u, self.np_u[1])

        if np_b_max < np_b_min:
            print('下缘偏心距不合适')
        elif np_u_max < np_u_min:
            print('上缘偏心距不合适')
        elif np_u_max < np_b_min or np_b_max < np_u_min:
            print('上下缘偏心距对应预应力没有交集！！！！')
            print(f'下缘最大/最小预应力(kN):{np_b_max:.0f} / {np_b_min:.0f}')
            print(f'下缘所需最多/最少根数:{np_b_max/pro_N_1:.0f} / {np_b_min/pro_N_1:.0f}')
            print(f'上缘最大/最小预应力(kN):{np_u_max:.0f} / {np_u_min:.0f}')
            print(f'上缘所需最多/最少根数:{np_u_max/pro_N_1:.0f} / {np_u_min/pro_N_1:.0f}')
        else:
            print(f'下缘最大/最小预应力(kN):{np_b_max:.0f} / {np_b_min:.0f}')
            print(f'下缘所需最多/最少根数:{np_b_max/pro_N_1:.0f} / {np_b_min/pro_N_1:.0f}')
            print(f'上缘最大/最小预应力(kN):{np_u_max:.0f} / {np_u_min:.0f}')
            print(f'上缘所需最多/最少根数:{np_u_max/pro_N_1:.0f} / {np_u_min/pro_N_1:.0f}')




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





