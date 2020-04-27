"""
此模块用于生成 mtc 命令流常用操作
适用 Midas 版本：15版本
"""

import pandas as pd

# 一些默认值 ====================================================================================================

units_default = 'KN, m, KJ, C'

# material_default = """
# 1, CONC, C50, 0, 0, , C, NO, 0.05, 1, JTG04(RC), , C50
# 2, STEEL, Strand1860, 0, 0, , C, NO, 0.02, 1, JTG04(S), , Strand1860
# 3, STEEL, Q345, 0, 0, , C, NO, 0.02, 1, GB12(S), , Q345
# """
material_default = """
1, CONC, C50, 0, 0, , C, NO, 0.05, 1, JTG04(RC), , C50
2, STEEL, Strand1860, 0, 0, , C, NO, 0.02, 1, JTG04(S), , Strand1860
3, CONC, C40, 0, 0, , C, NO, 0.05, 1, JTG04(RC), , C40
"""

stldcase_default = ('DEAD', 'TENDON', 'SUPERIMPOSED_DEAD', 'TEM_U', 'TEM_D')
loadgroup_default = ('dead', 'tendon', 'superimposed_deam', 'tem_u', 'tem_d')
bandgroup_default = ('support', 'rigid_connection', 'constraint', 'elastic_connection')
tendongroup_default = ('T1', 'T2')

tdm_defaul = (
    'NAME=C50, JTG, 50000, 67, 1, 5, 3',
    'NAME=C40, JTG, 40000, 67, 1, 5, 3'
)

tdm_link_default = ('1, C50,', '3, C40,')

num_default = [0, 0, 0, 0]
mct_list_default = []


# 初始化 mct 文件 ====================================================================================================

def usual_data(
        mct_list=mct_list_default,
        unit=units_default,
        material=material_default,
        stldcase=stldcase_default,
        loadgroup=loadgroup_default,
        bandgroup=bandgroup_default,
        tendongroup=tendongroup_default
):
    """对单位、材料、工况、荷载组等进行初始化"""
    mct_list.append('*UNIT')
    mct_list.append(unit)
    mct_list.append('*MATERIAL')
    mct_list.append(material.strip())
    mct_list.append('*TDM-TYPE')
    mct_list.append('\n'.join(tdm_defaul))
    mct_list.append('*TDM-LINK')
    mct_list.append('\n'.join(tdm_link_default))
    mct_list.append('*STLDCASE')
    stldcase = [i + ', USER,' for i in stldcase]
    mct_list.append('\n'.join(stldcase))
    mct_list.append('*LOAD-GROUP')
    mct_list.append('\n'.join(loadgroup))
    mct_list.append('*BNDR-GROUP')
    mct_list.append('\n'.join(bandgroup))
    mct_list.append('*TENDON-GROUP')
    mct_list.append('\n'.join(tendongroup))


# 定义截面 ====================================================================================================

def psv_sec_start(mct_list=mct_list_default):
    mct_list.append('*SECT-PSCVALUE')


def sec_start(mct_list=mct_list_default):
    mct_list.append('*SECTION')


def add_my_sec(sec, mct_list=mct_list_default, num=num_default):
    """添加截面（截面号自动添加）"""
    sec_add_num = []
    for i in sec:
        num[0] += 1
        sec_add_num.append('SECT = ' + str(num[0]) + '\n'.join(i))
    mct_list.append('\n'.join(sec_add_num))


# def add_midas_sec(name, loc='CC', what='SB', geo=(2, 2), mct_list=mct_list_default, num=num_default):
#     """添加 midas 常规截面"""
#     num[0] += 1
#     sec_str = f'{num[0]}, DBUSER, {name}, {loc}, 0, 0, 0, 0, 0, 0, YES,' \
#               f' {what}, 2, {geo[0]}, {geo[1]}, 0, 0, 0, 0, 0, 0, 0, 0'
#     mct_list.append(sec_str)


def add_change_sec(sec1, sec2, sec_all, name, pt='CT', mct_list=mct_list_default, num=num_default, hf=1, vf=1):
    """添加变截面"""
    num[0] += 1
    sec_1 = sec_all[sec1 - 1]
    sec_2 = sec_all[sec2 - 1]
    sec_1 = [i.replace('POLY=', 'POLY=YES,') for i in sec_1]
    sec_2 = [i.replace('POLY=', 'POLY=NO,') for i in sec_2]
    sec_str = [f'SECT={num[0]}, TAPERED, {name}, {pt}, 0, 0, 0, 0, 0, 0, 0, 0, YES, GEN, {hf}, {vf}, YES']
    sec_str += sec_1[1: 4]
    sec_str += sec_2[1: 4]
    sec_str += sec_1[4:]
    sec_str += sec_2[4:]
    mct_list.append('\n'.join(sec_str))
    return num[0]


def start_change_group(mct_list):
    """开始变截面组"""
    mct_list.append('*TS-GROUP')


def add_change_group(name, elem, z_change='', z_side='', z_d='', y_change='',
                     y_side='', y_d='', mct_list=mct_list_default):
    """添加变截面组"""
    z_type = 'LINEAR' if z_change == "" else 'QUADRATIC'
    y_type = 'LINEAR' if y_change == "" else 'QUADRATIC'

    mct_list.append(f'{name}, {elem}, {z_type}, {z_change}, {z_side}, {z_d}, '
                    f'{y_type}, {y_change}, {y_side}, {y_d}, 0')


# 节点生成单元函数 ====================================================================================================

def elem_from_nodes(nodes, s=1, what='nothing', m=1, b=0, num=num_default):
    """连续节点生成单元"""
    elem = pd.DataFrame(columns=['m', 's', 'b', 'n1', 'n2', 'what'])
    for i, j in enumerate(nodes[:-1]):
        num[2] += 1
        n1 = j
        n2 = nodes[i + 1]
        elem.loc[num[2]] = [m, s, b, n1, n2, what]
    return elem


def elem_form_double_nodes(double_nodes, s=1, what='nothing', m=1, b=0, num=num_default):
    """两两节点生成单元"""
    elem = pd.DataFrame(columns=['m', 's', 'b', 'n1', 'n2', 'what'])
    for i, j in zip(double_nodes[0], double_nodes[1]):
        n1 = i
        n2 = j
        num[2] += 1
        elem.loc[num[2]] = [m, s, b, n1, n2, what]
    return elem


def node_to_mct(node, mct_list=mct_list_default):
    """节点写入 mct"""
    mct_list.append('*NODE')
    for i, j in node.iterrows():
        mct_list.append(f'{i}, {j["x"]}, {j["y"]}, {j["z"]}')


def elem_to_mct(elem, elem_type='BEAM', mct_list=mct_list_default):
    """单元写入 mct"""
    mct_list.append('*ELEMENT')
    for i, j in elem.iterrows():
        mct_list.append(f"{i}, {elem_type}, {j['m']}, {j['s']}, {j['n1']}, {j['n2']}, {j['b']}")


def add_group(group_list, node_list, elem_list, mct_list=mct_list_default):
    """添加结构组"""
    mct_list.append('*GROUP')
    for i, j in enumerate(group_list):
        mct_list.append(f'{j}, {node_list[i]}, {elem_list[i]}, 0')


# 添加刚臂、边界条件等 ================================================================================================

def start_boundary(mct_list=mct_list_default):
    mct_list.append('*CONSTRAINT')


def start_link(mct_list=mct_list_default):
    mct_list.append('*ELASTICLINK')


def add_boundary(node, group, band_type, mct_list=mct_list_default):
    """
    对单个节点添加约束
    :param node: (int)，单元号
    :param group: (str)，约束组名称
    :param band_type: (str)，约束种类
    :param mct_list: (list)，mct列表
    :return:
    """
    mct_list.append(f'{node}, {band_type}, {group}')


def add_elastic_link(n1, n2, group, num=num_default, mct_list=mct_list_default,
                   sdx=0, sdy=0, sdz=0, srx=0, sry=0, srz=0):
    """
    添加弹性连接
    :param n1: 左节点
    :param n2: 右节点
    :param group: 边界组
    :param num: 编号
    :param mct_list: mct 列表
    :param sdx: x 线刚度（kN/m)
    :param sdy: y 线刚度
    :param sdz: z 线刚度
    :param srx: x 转角刚度
    :param sry: y 转角刚度
    :param srz: z 转角刚度
    :return:
    """
    num[3] += 1
    mct_list.append(f'{num[3]}, {n1}, {n2}, GEN, 0, '
                    f'{sdx}, {sdy}, {sdz}, {srx}, {sry}, {srz}, NO, 0.5, 0.5, {group}')


def add_rigid_link(n1, n2, group, num=num_default, mct_list=mct_list_default):
    """
    添加刚性连接
    :param n1: 左节点
    :param n2: 右节点
    :param group: 边界组
    :param num: 编号
    :param mct_list: mct 列表
    :return:
    """
    num[3] += 1
    mct_list.append(f'{num[3]}, {n1}, {n2}, RIGID,0, NO, 0.5, 0.5, {group}')


# 荷载 ================================================================================================

def start_stld(stld_name='', mct_list=mct_list_default):
    """开始某个工况"""
    mct_list.append(f'*USE-STLD, {stld_name}')


def add_self_weight(weight_factor=1, mct_list=mct_list_default, group=''):
    """
    添加自重荷载
    :param weight_factor: 自重系数，默认为 1
    :param mct_list: mct 列表
    :param group: 荷载组
    :return:
    """
    mct_list.append(f'*SELFWEIGHT, 0, 0, {-weight_factor}, {group}')


def start_node_load(mct_list=mct_list_default):
    """
    开始节点荷载
    :param mct_list:
    :return:
    """
    mct_list.append('*CONLOAD')


def node_load(node, mct_list=mct_list_default, group='', fx=0, fy=0, fz=0, mx=0, my=0, mz=0):
    """
    节点荷载
    :param node:
    :param mct_list:
    :param group:
    :param fx:
    :param fy:
    :param fz:
    :param mx:
    :param my:
    :param mz:
    :return:
    """
    mct_list.append(f'{node}, {fx}, {fy}, {fz}, {mx}, {my}, {mz}, {group}')


def start_tem_load(mct_list=mct_list_default):
    """
    开始温度荷载
    :param mct_list:
    :return:
    """
    mct_list.append('*ELTEMPER')


def tem_load(elem, tem, mct_list=mct_list_default, group=''):
    """
    添加温度荷载
    :param elem:
    :param tem:
    :param mct_list:
    :param group:
    :return:
    """
    mct_list.append(f'{elem}, {tem}, {group}')


def tendon_prop(name, tendon_a, material, mct_list=mct_list_default):
    """
    批量添加钢束特性
    :param name: 钢束名称（list）
    :param tendon_a: 钢束面积（list），单位 m2
    :param material: 材料号
    :param mct_list: mct 列表
    :return:
    """
    mct_list.append('*TDN-PROPERTY')
    for i, j in zip(name, tendon_a):
        mct_list.append(f'{i}, INTERNAL, {material}, {j}, 0.1, 2, 0, 0.3,'
                        f' 0.0066, 1.86326e+006, 1.56906e+006, POST, 0.006, 0.006, YES, 0, NO, 1, 1860000, 0, YES')


def start_tendon_type(mct_list=mct_list_default):
    mct_list.append('*TDN-PROFILE')


def tendon_type(
        type_name,
        prop_name,
        elem,
        tendon_n,
        origin_coordinate,
        x_direction,
        y_pt,
        z_pt,
        mct_list=mct_list_default
):
    """
    定义钢束线型
    :param type_name: 线型名称
    :param prop_name: 钢束属性名称
    :param elem: 单元列表（list）
    :param tendon_n: 钢束数量
    :param origin_coordinate: 原点坐标（x,y,z)
    :param x_direction: 局部坐标系 x 轴方向
    :param y_pt: y方向布置 [(x,y,r)...]
    :param z_pt: z 方向布置 [(x,z,r)...]
    :param mct_list: mct 列表
    :return:
    """
    elem = [str(i) for i in elem]
    mct_list.append(f'NAME={type_name}, {prop_name}, {" ".join(elem)}, 0, 0, ROUND, 2D')
    mct_list.append(f', AUTO1, , , YES, {tendon_n}')
    mct_list.append(f'STRAIGHT, {origin_coordinate[0]}, {origin_coordinate[1]},'
                    f'{origin_coordinate[2]}, {x_direction}, 0, 0')
    mct_list.append(f'0, YES, Y, 0')
    for i in y_pt:
        mct_list.append(f'Y={i[0]}, {i[1]}, NO, 0, {i[2]}, NONE, , , , ')
    for i in z_pt:
        mct_list.append(f'Z={i[0]}, {i[1]}, NO, 0, {i[2]}, NONE, , , , ')


def start_tendon_f(mct_list=mct_list_default):
    """开启钢束加载"""
    mct_list.append('*TDN-PRESTRESS')


def tendon_f(name, stress=1.395e6, mct_list=mct_list_default):
    """钢束加载"""
    mct_list.append(f'{name}, STRESS, BOTH, {stress}, {stress}, 0, tendon')


def start_com(mct_list=mct_list_default):
    """开始荷载组合"""
    mct_list.append('*LOADCOMB')


def com_f(name, stld_list, factor_list, type_list, type=0, mct_list=mct_list_default):
    """荷载组合"""
    mct_list.append(f'NAME={name}, GEN, ACTIVE, 0, {type}, , 0, 0')
    com_str = []
    for i, j, k in zip(stld_list, factor_list, type_list):
        com_str.append(f'{k}, {i}, {j}')
    com_str = ','.join(com_str)
    mct_list.append(com_str)


def start_stage(mct_list=mct_list_default):
    """开始施工阶段"""
    mct_list.append('*STAGE')


def add_stage(name, time, elem=[], elem_time=[], banr=[], load=[], mct_list=mct_list_default):
    mct_list.append(f'NAME={name}, {time}, YES, NO')
    if elem:
        elem_str = []
        for i, j in zip(elem, elem_time):
            elem_str.append(f'{i}, {j}')
        mct_list.append('AELEM=' + ','.join(elem_str))
    if banr:
        banr_str = []
        for i in banr:
            banr_str.append(f'{i}, ORIGINAL')
        mct_list.append(f'ABNDR={",".join(banr_str)}')
    if load:
        load_str = []
        for i in load:
            load_str.append(f'{i}, FIRST')
        mct_list.append('ALOAD=' + ','.join(load_str))


