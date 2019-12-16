import math
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

Points = ['0']
Flat = []
point_coordinate = []
points_draw = np.zeros([6, 3])


def make_point(x, y, z):            #根据X、Y、Z坐标点生成向量
    return np.array([x, y, z])


def volume(a, b, c, d):             #返回四面体体积VABCD
    return math.fabs(np.dot(np.cross((b - a), (c - a)), (d - a))) / 6.000


def distance(a, b):
    dist = np.linalg.norm(a - b)
    return dist


def find_center(a, b, c, d, e):
    global G
    V_ABCD = volume(a, b, c, d)
    V_ABCE = volume(a, b, c, e)

    G_ABCD = (a + b + c + d) / 4.000
    G_ABCE = (a + b + c + e) / 4.000

    G = G_ABCD + ((G_ABCE - G_ABCD) * V_ABCE) / (V_ABCD + V_ABCE)


def same_side(p1, p2, p3):
    global Flat
    normal_vector = np.cross((p2 - p1), (p3 - p1))
    plus = False
    minus = False
    for i in range(1, 6):
        if (Points[i] == p1).all() or (Points[i] == p2).all() or (Points[i] == p3).all():
            continue
        temp = np.dot((Points[i] - p1), normal_vector)
        if temp < 0:
            minus = True
        elif temp > 0:
            plus = True
    if plus and minus:
        return False
    Flat.clear()
    Flat.append(p1)
    Flat.append(p2)
    Flat.append(p3)
    for i in range(1, 6):
        if (Points[i] == p1).all() or (Points[i] == p2).all() or (Points[i] == p3).all() :
            continue
        temp = np.dot((Points[i] - p1), normal_vector)
        if temp == 0:
            Flat.append(Points[i])
    return True


def shadow(p, p1, p2, p3):
    normal_vector = np.cross((p2 - p1), (p3 - p1))
    normal_vector0 = np.dot(normal_vector, np.dot(normal_vector, (p - p1))) / np.linalg.norm(normal_vector) / np.linalg.norm(normal_vector)
    return p - normal_vector0


def shadow_height(p, p1, p2, p3):
    normal_vector = np.cross((p2 - p1), (p3 - p1))
    return math.fabs(np.dot(normal_vector, (p - p1)) / np.linalg.norm(normal_vector))


def in_the_middle(a, b, c, p):
    ab = b - a
    ac = c - a
    ap = p - a
    v1 = np.cross(ab, ac)
    v2 = np.cross(ab, ap)
    return np.dot(v1, v2) >= 0


def point_in_triangle(a, b, c, p):
    return in_the_middle(a, b, c, p) and in_the_middle(b, c, a, p) and in_the_middle(c, a, b, p)


def point_in_area(p):
    for i in range(0, len(Flat)):
        for j in range(i + 1, len(Flat)):
            for k in range(j + 1, len(Flat)):
                if point_in_triangle(Flat[i], Flat[j], Flat[k], p):
                    return True
                return False


def distance_to_line(p, a, b):
    if np.dot((b - a), (p - a)) > 0.000 and np.dot((a - b), (p - b)) > 0.000:
        return math.sqrt(distance(a, p) ** 2 - (np.dot((b - a), (p - a)) / np.linalg.norm(b - a)) ** 2)
    return min(distance(a, p), distance(b, p))


def balance(p):
    if not point_in_area(p):
        return False
    for i in range(0, len(Flat)):
        for j in range(i + 1, len(Flat)):
            p_i = Flat[i]
            p_j = Flat[j]
            for k in range(0, len(Flat)):
                if k != i and k != j:
                    break
            p_k = Flat[k]

            plus = False
            minus = False
            normal_vector = np.cross(np.cross((p_j - p_i), (p_k - p_i)),(p_j - p_i))
            for l in range(0, len(Flat)):
                if l != i and l != j:
                    if np.dot(normal_vector, (Flat[l] - p_i)) > 0:
                        plus = True
                    elif np.dot(normal_vector, (Flat[l] - p_i)) < 0:
                        minus = True
            if plus and minus:
                continue
            if distance_to_line(p, p_i, p_j) < 0.2:
                return False
        return True


def find_hull():
    min_answer = 2000.000
    max_answer = 0.000
    for i in range(1, 6):
        for j in range(i + 1, 6):
            for k in range(j + 1, 6):
                if same_side(Points[i], Points[j], Points[k]):
                    G0 = shadow(G, Points[i], Points[j], Points[k])
                    if balance(G0):
                        min_answer = min(min_answer, shadow_height(F, Points[i], Points[j], Points[k]) )
                        max_answer = max(max_answer, shadow_height(F, Points[i], Points[j], Points[k]))
    print('%.5f' % min_answer + " " + '%.5f' % max_answer + "\n")


def input_pts(pts):
    # pts = input("Input Points >>")
    pts = np.array(point_coordinate).reshape(6, 3)

    # 输入小于1000

    # ABC不共线

    # DE分布ABC两侧

    # F必须严格在两四面体内部

    return pts


def output(pts):
    # 用123456代表6个点，输出二维矩阵[[点1，点2，点3，距离]，...]
    return


def draw(pts):
    x = pts[:, 0]
    y = pts[:, 1]
    z = pts[:, 2]
    label = np.array(['A', 'B', 'C', 'D', 'E', 'F'])
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x[:-1], y[:-1], z[:-1], color='b')
    ax.scatter(x[-1], y[-1], z[-1], color='r')

    for i in range(len(x)):
        ax.text(x[i], y[i], z[i], label[i])

    for i in range(len(x) - 1):
        for j in range(len(x) - 1):
            if i not in [3, 4] or j not in [3, 4]:
                ax.plot([x[i], x[j]], [y[i], y[j]], [z[i], z[j]], c='g')

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.view_init(elev=20, azim=35)
    # ax.grid(False)
    plt.show()


def main():
    global A, B, C, D, E, F, G
    global points_draw,point_coordinate
    case = 0
    while True:
        string = input("请输入坐标，以空格隔开：")
        if string == "0":
            break
        point_coordinate = string.split(" ")
        point_coordinate = [int(x) for x in point_coordinate]
        # ------
        points_draw = input_pts(points_draw)
        draw(points_draw)
        print(points_draw)
        # --------
        point_coordinate.append(0)
        point_coordinate.append(0)
        point_coordinate.append(0)
        creatVar = globals()
        for name_ascii in range(65, 72):
            index = name_ascii - 65
            creatVar[chr(name_ascii)] = make_point(point_coordinate[index * 3], point_coordinate[index * 3 + 1],
                                                   point_coordinate[index * 3 + 2])
            Points.append(eval(chr(name_ascii)))
        print("Case " + str(case) + ":")
        find_center(A, B, C, D, E)
        find_hull()


if __name__ == '__main__':
    main()
