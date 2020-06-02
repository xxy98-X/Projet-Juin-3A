import numpy

class Element(object):
    def __init__(self, node_label, A, E, L, angle=None, cosine=None):
        self.node_label = node_label
        self.A = A
        self.E = E
        self.L = L
        self.k = self.A * self.E / self.L

        if angle != None and cosine == None:
            C, S, CS = self.__calculate_CS(angle)
        elif cosine != None and angle == None:
            C = cosine
            S = numpy.sqrt(1-C**2)
            CS = C * S
        self.C, self.S = C, S
        self.kk = self.__generate_local_matrix(C,S,CS,self.k)

    def __generate_local_matrix(self, C, S, CS, K):
        kk = numpy.array([[0 for i in range(4)] for i in range(4)])
        kk[0][0] = C ** 2 * K
        kk[2][2] = C ** 2 * K
        kk[1][1] = S ** 2 * K
        kk[3][3] = S ** 2 * K
        kk[0][1] = CS * K
        kk[1][0] = CS * K
        kk[2][3] = CS * K
        kk[3][2] = CS * K
        kk[0][2] = - C ** 2 * K
        kk[2][0] = - C ** 2 * K
        kk[1][3] = - S ** 2 * K
        kk[3][1] = - S ** 2 * K
        kk[0][3] = - CS * K
        kk[1][2] = - CS * K
        kk[2][1] = - CS * K
        kk[3][0] = - CS * K
        return kk

    def __calculate_CS(self, angle):
        C = numpy.cos(angle)
        S = numpy.sin(angle)
        CS = C*S
        return C, S, CS






def generate_global_matrix(K, ElementSet):
    for e in ElementSet:
        i, j =e.node_label
        kk = e.kk
        row = (i - 1) * 2
        col = (j - 1) * 2
        for m in range(kk.shape[0]):
            for n in range(kk.shape[1]):
                if m < 2:
                    if n < 2:
                        K[row + m % 2][row + n % 2] += kk[m][n]
                        continue
                    else:
                        K[col + m % 2][row + n % 2] += kk[m][n]
                        continue
                else:
                    if n < 2:
                        K[row + m % 2][col + n % 2] += kk[m][n]
                        continue
                    else:
                        K[col + m % 2][col + n % 2] += kk[m][n]
                        continue

    return K

def calculate_node_displacement(KK, U, F):
    index = numpy.where(U==0)
    new_U = numpy.delete(U, index, 0)
    KK_del_row = numpy.delete(KK, index, 0)
    KK_del_col = numpy.delete(KK_del_row, index, 1)
    print("kk_del_col={}".format(KK_del_col))
    new_F = numpy.delete(F, index, 0)
    new_KK = numpy.array([[-1, KK_del_col[0][1]],[0,KK_del_col[1][1]]])
    new_f = numpy.array([-new_U[0]*KK_del_col[0][0], new_F[1]-new_U[0]*KK_del_col[1][0]])

    u_f, u_u = numpy.linalg.solve(new_KK, new_f)


    U[numpy.where(U == "u_unknown")] = u_u
    print("u_f={}".format(u_f))
    F = KK.dot(U)
    return U,F



if __name__ == '__main__':
    ElementSet = []
    E1 = Element((1, 2), 6e-4, 210e6, 5, cosine=0.6)
    print("E1.kk={}".format(E1.kk))
    E2 = Element((1, 3), 6e-4, 210e6, 4, angle=numpy.pi/2)
    print("E2.kk={}".format(E2.kk))
    ElementSet.append(E1)
    ElementSet.append(E2)
    K = numpy.array([[0 for i in range(6)] for i in range(6)])
    new_KK = generate_global_matrix(K, ElementSet)
    print("new_KK={}".format(new_KK))
    U = numpy.array([-0.05, "u_unknown", 0, 0, 0, 0], dtype=object)
    F = numpy.array(["f_unknown", 1000, "f_unknown", "f_unknown", "f_unknown", "f_unknown"], dtype=object)
    U, F = calculate_node_displacement(new_KK, U, F)
    print("U={}".format(U))
    print("F={}".format(F))

for e in ElementSet:
    node_i, node_j = e.node_label
    C, S = e.C, e.S
    print("C={0};S={1}".format(C, S))
    T = numpy.array([[C, S, 0, 0], [0, 0, C, S]])
    new_U = U[[(node_i-1)*2, (node_i-1)*2+1, (node_j-1)*2, (node_j-1)*2+1]]
    print("new_U={}".format(new_U))
    local_U = T.dot(new_U)
    local_F = e.k*numpy.array([[1, -1],[-1, 1]]).dot(local_U)
    print("local_F={}".format(local_F))


