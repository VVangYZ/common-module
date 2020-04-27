"""
用于处理 cad 图形
1. 截面计算（单位：m）
"""

import win32com.client as win32
import numpy as np
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

# acad = win32.Dispatch("AutoCAD.Application")
# doc = acad.ActiveDocument
# ms = doc.ModelSpace


def find_pt(ps, p1=(0, 0), relation='min'):
    """
    用于确定一组点中距离某点最远（或最近）的点
    :param ps: 一组点
    :param p1: 定点，默认为 (0, 0)
    :param relation: 'max' or 'min' ('min' for default)
    :return: 距离最近或最远的点坐标
    """
    ps = np.array(ps)
    p1 = np.array(p1)
    ls = np.linalg.norm(p1 - ps, axis=1, keepdims=True).flatten()
    order = np.argsort(ls)
    if relation == 'min':
        return ps[order[0]]
    elif relation == 'max':
        return ps[order[-1]]
    else:
        print('wrong!')


class ClosedLine:
    """表示一根封闭多段线"""

    def __init__(self, pl):
        id_origin = pl.coordinates
        self.id = np.array(id_origin).reshape(len(id_origin) // 2, 2)
        self.area = pl.area
        self.length = pl.length
        self.layer = pl.layer


class OpenLine:
    """表示一根开放的多段线"""

    def __init__(self, pl):
        id_origin = pl.coordinates
        self.id = np.array(id_origin).reshape(len(id_origin) // 2, 2)
        self.length = pl.length
        self.layer = pl.layer
        self.area = 0


class OneSec:
    """表示一个截面"""

    def __init__(self, pls):
        """
        单独截面
        :param pls: 多条 cad 多段线对象组成的列表，其中应有一根指示控制点的线
        """

        # 原始线和控制点
        self.pls_origin = []
        self.points = []
        for i in pls:
            if i.area == 0:
                self.points = i.id
            else:
                self.pls_origin.append(i)

        # 对线段进行排序，得到分离节点和完整节点
        self.areas = np.array([i.area for i in self.pls_origin])
        self.pls = np.array(self.pls_origin)[np.argsort(self.areas)][::-1]
        self.ids_sep = [i.id for i in self.pls]
        self.ids = [j.tolist() for i in self.ids_sep for j in i]

        # 获取节点连接方式
        self.faces = []
        id_num = 0
        for i in self.ids_sep:
            id_num_0 = id_num
            for j in i:
                connect = [id_num, id_num_0] if id_num + 1 == id_num_0 + len(i) else [id_num, id_num + 1]
                self.faces.append(connect)
                id_num += 1

        # 定义其他所需值
        self.geo = 0
        self.mesh = 0
        self.sec = 0
        self.prop = {}
        self.stress = 0
        self.corner = []
        self.ids_to_c = []

    def sec_cal(self, mesh=0.01, d=0.03):
        """
        对单个截面进行属性计算
        :param mesh: 截面划分单元尺寸
        :param d: 截取形心附近应力范围
        :return: 无
        """

        self.geo = sections.CustomSection(self.ids, self.faces, self.points[1:], [self.points[0]])
        self.mesh = self.geo.create_mesh(mesh_sizes=[mesh])
        self.sec = CrossSection(self.geo, self.mesh)
        self.sec.plot_mesh()
        self.sec.calculate_geometric_properties()
        self.sec.calculate_warping_properties()

        # 获取截面属性
        prop = self.sec.section_props
        self.prop['center'] = self.sec.get_c()
        self.ids_to_c = [i - self.prop['center'] for i in self.ids_sep]
        self.prop['area'] = prop.area
        self.prop['as'] = [prop.A_s22, prop.A_s11]
        self.prop['i'] = [prop.j, prop.ixx_c, prop.iyy_c]
        pts = np.array(self.ids)
        left = prop.cx - pts[:, 0].min()
        right = pts[:, 0].max() - prop.cx
        top = pts[:, 1].max() - prop.cy
        bot = prop.cy - pts[:, 1].min()
        self.prop['c'] = [right, left, top, bot]

        self.stress = self.sec.calculate_stress(Vx=1, Vy=1)
        stresses = self.stress.get_stress()
        dy = self.sec.get_c()[1] - self.sec.mesh_nodes[:, 1]
        dx = self.sec.get_c()[0] - self.sec.mesh_nodes[:, 0]
        qyb = stresses[0]['sig_zy_vy'][dx < d].max() * prop.ixx_c
        qzb = stresses[0]['sig_zx_vx'][dy < d].max() * prop.iyy_c
        self.prop['q'] = [qyb, qzb]
        self.prop['p'] = [self.pls[0].length, sum([i.length for i in self.pls[1:]])]

        # 获取角点
        pt_all = self.ids_to_c[0]
        pt_1 = pt_all[(pt_all[:, 0] < 0) & (pt_all[:, 1] > 0)]
        pt_2 = pt_all[(pt_all[:, 0] > 0) & (pt_all[:, 1] > 0)]
        pt_3 = pt_all[(pt_all[:, 0] < 0) & (pt_all[:, 1] < 0)]
        pt_4 = pt_all[(pt_all[:, 0] > 0) & (pt_all[:, 1] < 0)]
        pt_1 = find_pt(pt_1, relation='max')
        pt_2 = find_pt(pt_2, relation='max')
        pt_3 = find_pt(pt_3, relation='max')
        pt_4 = find_pt(pt_4, relation='max')
        self.corner = [pt_1, pt_2, pt_4, pt_3]


class CadSecs:
    """表示一组 cad 截面"""

    def __init__(self, cad_ms):
        """
        定义一组 cad 截面的类
        :param cad_ms: cad 文件的 model space，包含多个图层，每个图层中绘制截面
        """
        # 实例化各多段线
        self.cls = []
        self.ols = []
        for i in range(cad_ms.count):
            the_item = cad_ms.item(i)
            if the_item.objectname == 'AcDbPolyline':
                if the_item.closed:
                    self.cls.append(ClosedLine(the_item))
                else:
                    self.ols.append(OpenLine(the_item))

        # 多段线分组
        self.sec_name = np.unique([i.layer for i in self.cls])
        self.sec_lines = {}
        for i in self.sec_name:
            self.sec_lines.setdefault(i, [])
            for j in self.cls:
                if j.layer == i:
                    self.sec_lines[i].append(j)
            for j in self.ols:
                if j.layer == i:
                    self.sec_lines[i].append(j)

        # 实例化各截面
        self.secs = []
        for i in self.sec_lines:
            self.secs.append(OneSec(self.sec_lines[i]))
            self.secs[-1].sec_cal()

    def secs_to_midas(self, pt='CT'):
        """将各个截面数据转化为 midas 命令流"""
        midas_str = []
        for m, i in enumerate(self.secs):
            stiff1 = [i.prop['area']] + i.prop['as'] + i.prop['i']
            stiff2 = i.prop['c'] + i.prop['q'] + i.prop['p'] + [i.prop['c'][0], i.prop['c'][3]]
            stiff3 = [j[0] for j in i.corner] + [j[1] for j in i.corner]
            stiff = [stiff1, stiff2, stiff3]
            stiff_str = [f', VALUE, {self.sec_name[m]}, {pt}, 0, 0, 0, 0, 0, 0, YES, GEN, YES, YES']
            for j in stiff:
                stiff_str.append(','.join([str(round(k, 4)) for k in j]))

            for j, k in enumerate(sorted(i.ids_to_c, key=lambda x: x[:, 0].min())):
                # poly = k.flatten()
                # poly_str = ','.join([str(round(m, 4)) for m in poly])
                if j == 0:
                    poly = k[::-1].flatten()
                    poly_str = ','.join([str(round(m, 4)) for m in poly])
                    stiff_str.append('OPOLY=' + poly_str)
                # elif j == len(i.ids_to_c) - 1:
                #     poly = k[::-1].flatten()
                #     poly_str = ','.join([str(round(m, 4)) for m in poly])
                #     stiff_str.append('IPOLY=' + poly_str)
                else:
                    poly = k.flatten()
                    poly_str = ','.join([str(round(m, 4)) for m in poly])
                    stiff_str.append('IPOLY=' + poly_str)


            midas_str.append(stiff_str)

        return midas_str


# a = CadSecs(ms)
# b = a.secs_to_midas()
# a.secs[0].sec.plot_mesh()
# print(b[0])