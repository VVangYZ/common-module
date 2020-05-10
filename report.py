"""用于根据计算结果生成计算报告

包含：风荷载计算
"""

from pylatex import Document, Section, Subsection, Math, Command, Package, Enumerate, Tabular, LongTabu, Figure, LongTable
from pylatex.utils import NoEscape, bold
import pandas as pd
import pile as pl
import numpy as np

geometry_options = 'a4paper, left=2.5cm, right=2.5cm, top=3cm, bottom=3cm'
document_options = 'cs4size, UTF8'

class Report():
    """生成报告基本模板
    
    参数：
        name：文档名称
        geo：文档页面设置，默认为 geometry_options
        doo：文档选项，默认为 document_options，小四字体为正文字体
    """

    def __init__(self, name, geo=geometry_options, doo=document_options):
        self.document_options = doo
        self.geometry_options = geo
        self.doc = Document(name, documentclass='ctexart',document_options=self.document_options, geometry_options=self.geometry_options)
        # 加入宏包
        self.doc.packages.append(Package("xeCJK"))
        self.doc.packages.append(Package("fontspec"))
    
    def set_format1(self):
        """设置为常用格式，英文字体 Times，章节标题左对齐，大标题"""

        self.doc.preamble.append(Command('setmainfont', 'Times New Roman'))
        self.doc.preamble.append(Command('CTEXsetup', 'section', NoEscape(r'format={\raggedright\bfseries\Large}')))

    def set_format2(self):
        """设置为常用格式，英文字体 Times，章节标题左对齐，小标题"""

        # self.doc.preamble.append(Command('setmainfont', 'Times New Roman'))
        self.doc.preamble.append(Command('CTEXsetup', 'section', NoEscape(r'format={\raggedright\bfseries\large}')))
        self.doc.preamble.append(Command('pagestyle', 'plain'))
        self.doc.packages.append(Package("mathptmx"))   # Times 字体宏包
    
    def tex_out(self):
        '''得到计算报告'''

        self.doc.generate_tex()
    

class ReportWind(Report):
    """桥梁静风荷载计算报告

    """

    def __init__(self, name, geo=geometry_options, doo=document_options):
        super().__init__(name, geo, doo)
        self.doc.preamble.append(Command('title', NoEscape(r'\heiti 桥梁风荷载计算')))
        self.doc.preamble.append(Command('author', 'WYZ'))
        self.doc.preamble.append(Command('date', NoEscape(r'\today')))
        self.doc.append(NoEscape(r'\maketitle'))
    
    def wind_report(self, wind_beam, wind_pier):
        """添加风荷载计算内容"""

        self.wind = wind_beam.wind

        t1 = f'''桥梁抗风风险区域：R{self.wind.R}\n
        桥位地表分类：{self.wind.surface_class}\n
        十年重现期风作用水平：W1={self.wind.W1:.3f} m/s\n
        百年重现期风作用水平：W2={self.wind.W2:.3f} m/s
        '''
        t2 = '根据规范4.1可得，基本风速值为'
        m1 = [f'U_{{10}}={self.wind.W2:.3f}', '\\,m/s']

        t3 = '根据规范4.2.1可得，地表相关参数为'
        m2 = [f'\\alpha_0={self.wind.alpha0:.2f}', '\\quad', f'z_0={self.wind.z0:.2f}']

        t5 = '根据规范4.2.4可得，桥梁设计基本风速为'
        m4 = ['U_{s10}=k_cU_{10}=', f'{self.wind.kc:.3f}\\times{self.wind.U10:.3f}={self.wind.Us10:.3f}', '\\,m/s']


        t4 = '根据规范4.2.2可得，主梁基准高度为'
        m3 = [f'Z={wind_beam.z:.2f}', '\\,m']

        t6 = '根据规范4.2.6可得，主梁构件基准高度处的设计基准风速为'
        m5 = ['U_d=k_f\\left(\\frac{Z}{10}\\right)^{\\alpha_0}U_{s10}=',
        f'{self.wind.kf:.2f}\\times\\left(\\frac{{{wind_beam.z:.2f}}}{{10}}\\right)^{{{self.wind.alpha0:.2f}}}\\times{self.wind.Us10:.3f}=',
        f'{wind_beam.Ud:.3f}',
        '\\,m/s']

        # t7 = '根据规范4.2.9可得，施工阶段的设计风速为'
        # m6 = ['U_{sd}=k_{sf}U_{d}=', f'{self.wind.ksf}']

        t8 = '根据规范5.2.1可得，等效静阵风风速为'
        m7 = ['U_g=G_VU_d=', f'{wind_beam.GV:.2f}\\times{wind_beam.Ud:.3f}={wind_beam.Ug:.3f}', '\\,m/s']

        t9 = '根据规范5.3.1可得，横桥向风作用下主梁单位长度上的顺风向等效静阵风荷载为'
        m8 = ['F_g=\\frac{1}{2}\\rho U_g^2C_HD=',
        f'0.5\\times{wind_beam.rho:.2f}\\times{wind_beam.Ug:.3f}^2\\times{wind_beam.CH:.3f}\\times{wind_beam.D:.3f}=',
        f'{wind_beam.Fg:.3f}',
        '\\,N/m']


        with self.doc.create(Section('地质基本情况')):
            self.doc.append(NoEscape(t1))
        with self.doc.create(Section('风速参数')):
            self.doc.append(t2)
            self.doc.append(Math(data=m1, escape=False))
            self.doc.append(NoEscape('\n'))
            self.doc.append(t3)
            self.doc.append(Math(data=m2, escape=False))
            self.doc.append(NoEscape('\n'))
            self.doc.append(t5)
            self.doc.append(Math(data=m4, escape=False))
        with self.doc.create(Section('主梁风荷载')):
            self.doc.append(t4)
            self.doc.append(Math(data=m3, escape=False))
            self.doc.append(NoEscape('\n'))
            self.doc.append(t6)
            self.doc.append(Math(data=m5, escape=False))
            self.doc.append(NoEscape('\n'))
            self.doc.append(t8)
            self.doc.append(Math(data=m7, escape=False))
            self.doc.append(NoEscape('\n'))
            self.doc.append(t9)
            self.doc.append(Math(data=m8, escape=False))


class ReportPile(Report):
    """桥梁桩长计算报告

    """

    def __init__(self, name, geo=geometry_options, doo=document_options):
        super().__init__(name, geo, doo)
        self.doc.preamble.append(Command('title', NoEscape(r'\heiti 桥梁桩长计算')))
        self.doc.preamble.append(Command('author', ''))
        self.doc.preamble.append(Command('date', NoEscape(r'')))
        self.doc.append(NoEscape(r'\maketitle'))
    
    def pile_report(self, pile, pile2, F, factor=1.25):
        """添加桩基长度计算内容"""

        self.soil = pile.soil
        llist = pl.get_l(self.soil, pile.d, F, factor, pile.h1, pile.rho, pile.t, pile.k2, pile2.type, pile2.completion)
        ra_1, ra_2 = pl.pile_l(llist, self.soil, pile.d, F, pile.h1, pile.rho, pile.t, pile.k2, pile2.type, pile2.completion)

        if ra_1.max() > ra_2.max():
            ptype = '摩擦桩'
            ra = ra_1.max()
        else: 
            ptype = '端承桩'
            ra = ra_2.max()

        t1 = f'''桩基直径：$d={pile.d:.2f}\,m$\n
        桩基周长：$u={pile.u:.2f}\,m$\n
        桩基截面积：$A_p={pile.pd:.2f}\,m^2$\n
        桩基密度：$\gamma={pile.rho:.1f}\,kN/m^3$\n
        容许承载力随深度的修正系数：$k_2={pile.k2:.1f}$\n
        各土层加权平均重度：$\gamma_2={self.soil.rho:.1f}\,kN/m^3$\n
        清底系数：$m_0={pile.m0:.1f}$
        '''

        t2 = '根据规范5.3.3可得，摩擦桩单桩承载力为'
        m1 = ['[R_a]', '=\\frac{1}{2}u\\sum_{i=1}^nq_{ik}l_i+A_pq_r']
        t3 = '根据规范5.3.4可得，端承桩单桩承载力为'
        m2 = ['[R_a]=', 'c_1A_pf_{rk}', '+u\\sum_{i=1}^mc_{2i}h_if_{rki}', '+\\frac{1}{2}\\xi_su\sum_{i=1}^nl_iq_{ik}']
        t4 = '考虑桩身自重与置换土重，桩基承载力为'
        m3 = ['R_a', '=[R_a]-G_p+G_s']
        t5 = '代入不同长度桩长，可得摩擦桩与端承桩承载力如下图所示'
        t6 = '不同桩长具体承载力如下表所示'
        t7 = f'由上述分析可知，当桩长为{max(llist)}m时，{ptype}承载力为{ra:.0f}kN，安全系数为{ra/F:.2f}，桩基承载力可满足规范要求。'

        with self.doc.create(Section('地基基本情况')):
            with self.doc.create(LongTabu("p{4cm}XXXXX")) as soil_table:
                header_row1 = self.soil.prop.columns.to_list()[:-1]
                soil_table.add_hline()
                soil_table.add_row(header_row1, mapper=[bold])
                soil_table.add_hline()
                for i in self.soil.prop.index:
                    soil_table.add_row(self.soil.prop.iloc[i].to_list()[:-1])
                soil_table.add_hline()

        with self.doc.create(Section('桩基及其他参数取值')):
            self.doc.append(NoEscape(t1))
        
        with self.doc.create(Section('桩长计算')):
            self.doc.append(t2)
            self.doc.append(Math(data=m1, escape=False))
            self.doc.append(NoEscape('\n'))
            self.doc.append(t3)
            self.doc.append(Math(data=m2, escape=False))
            self.doc.append(NoEscape('\n'))
            self.doc.append(t4)
            self.doc.append(Math(data=m3, escape=False))
            self.doc.append(NoEscape('\n'))
            self.doc.append(t5)
            with self.doc.create(Figure(position='htbp')) as plot:
                plot.add_plot(width=NoEscape(r'1\textwidth'), dpi=300)
            self.doc.append(NoEscape('\n'))
            self.doc.append(t6)
            with self.doc.create(LongTable('p{1.5cm}|ll|ll')) as pll:
                pll.add_hline()
                pll.add_row(['桩长', '摩擦桩承载力', '安全系数', '端承桩承载力', '安全系数'])
                pll.add_hline()
                for i, j in enumerate(llist):
                    pll.add_row([j, f'{ra_1[i]:.0f}', f'{ra_1[i]/F:.2f}', f'{ra_2[i]:.0f}', f'{ra_2[i]/F:.2f}'])
                pll.add_hline()
            self.doc.append(t7)
            






            


