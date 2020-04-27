import clr
import numpy as np
import pandas as pd
import pymysql

clr.FindAssembly("SmartRoadBridge.Alignment")
clr.AddReference('SmartRoadBridge.Alignment')
clr.AddReference('System.Collections')

from SmartRoadBridge.Alignment import Align
from System import Array, String, Double
from System.Collections.Generic import List

# conn=pymysql.connect(host='cdb-2ashfo5g.bj.tencentcdb.com',
#     user = 'liu' ,passwd='liuning1234' ,port= 10033 ,db='testdb' ,charset='utf8')
# cur = conn.cursor()
# sql="select * from ei_tbl"
# cur.execute(sql)
# data = cur.fetchall()
# cur.close()
# conn.close()
# aa=np.array(data)

# AlignList={}

# for item in aa:
#     name=item[0].tolist()
#     icd = Array[String](item[1].split('\r\n'))
#     sqx = Array[String](item[2].split('\r\n'))
#     dmx = Array[String](item[3].split('\r\n'))
#     cg = Array[String](item[4].split('\r\n'))
#     AlignList[name]=Align(name,icd,sqx,dmx,cg)


def get_ei(ei_row, AlignList):
    name = ei_row.Name
    icd = Array[String](ei_row.ICD.strip().split('\r\n'))
    sqx = Array[String](ei_row.SQX.strip().split('\r\n'))
    dmx = Array[String](ei_row.DMX.strip().split('\r\n'))
    cg = Array[String](ei_row.CG.strip().split('\r\n'))
    AlignList[name] = Align(name,icd,sqx,dmx,cg)


def get_distance(rmp, sta, AlignList, ml='M1k'):
    """
    通过主线桩号，获取主线到匝道距离
    :param rmp: 匝道名称
    :param sta: 主线桩号
    :param ml: 主线名称，默认为 M1K
    :return:
    """
    rmp = rmp.upper()
    ml = ml.upper()

    cord1 = AlignList[ml].curPQX.GetCoord(sta)  # 主线坐标
    cord2 = AlignList[ml].curPQX.GetDir(sta)  # 主线切线坐标

    cord1 = [i for i in cord1]
    cord2 = [cord1[0] - cord2[1], cord1[1] + cord2[0]]  # 主线法线坐标
    # cord2 = [- cord2[1], cord1[0]]  # 主线法线坐标

    sta1 = AlignList[rmp].curPQX.GetStation(cord1[0], cord1[1], cord2[0], cord2[1])  # 匝道里程

    cord3 = AlignList[rmp].curPQX.GetCoord(sta1)  # 匝道切点
    cord3 = [i for i in cord3]
    dis1 = sum([(i[1] - i[0]) ** 2 for i in zip(cord1, cord3)]) ** 0.5  # 匝道切点和主线点距离

    return sta1, dis1


# print(get_distance('cca', 16430))

