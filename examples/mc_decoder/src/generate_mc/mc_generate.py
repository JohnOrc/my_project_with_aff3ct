# -*- coding: utf-8 -*-
# @Time    : 2020/12/17 9:35
# @Author  : zyl
# @File    : mc_generate.py
# @Contact : hit_zyl@126.com
from math import log2, ceil


def leqNum(l_list, L):
    ans = 0
    # 1 bit
    if len(l_list) >= 1:
        for i0 in range(L):
            if i0 <= l_list[-1]:
                ans = ans + 1
                if ans > L:
                    return ans
    # 2 bits
    if len(l_list) >= 2:
        for i0 in range(L):
            for i1 in range(i0 + 1, L):
                if i0 <= l_list[-2] and i1 <= l_list[-1]:
                    ans = ans + 1
                    if ans > L:
                        return ans
    # 3 bits
    if len(l_list) >= 3:
        for i0 in range(L):
            for i1 in range(i0 + 1, L):
                for i2 in range(i1 + 1, L):
                    if i0 <= l_list[-3] and i1 <= l_list[-2] and i2 <= l_list[-1]:
                        ans = ans + 1
                        if ans > L:
                            return ans
    # 4 bits
    if len(l_list) >= 4:
        for i0 in range(L):
            for i1 in range(i0 + 1, L):
                for i2 in range(i1 + 1, L):
                    for i3 in range(i2 + 1, L):
                        if i0 <= l_list[-4] and i1 <= l_list[-3] and i2 <= l_list[-2] and i3 <= l_list[-1]:
                            ans = ans + 1
                            if ans > L:
                                return ans
    # 5 bits
    if len(l_list) >= 5:
        for i0 in range(L):
            for i1 in range(i0 + 1, L):
                for i2 in range(i1 + 1, L):
                    for i3 in range(i2 + 1, L):
                        for i4 in range(i3 + 1, L):
                            if i0 <= l_list[-5] and i1 <= l_list[-4] and i2 <= l_list[-3] and i3 <= l_list[-2] and i4 <= l_list[-1]:
                                ans = ans + 1
                                if ans > L:
                                    return ans                            
    # 6 bits
    if len(l_list) >= 6:
        for i0 in range(L):
            for i1 in range(i0 + 1, L):
                for i2 in range(i1 + 1, L):
                    for i3 in range(i2 + 1, L):
                        for i4 in range(i3 + 1, L):
                            for i5 in range(i4 + 1, L):
                                if i0 <= l_list[-6] and i1 <= l_list[-5] and i2 <= l_list[-4] and i3 <= l_list[-3] and i4 <= \
                                        l_list[-2] and i5 <= l_list[-1]:
                                    ans = ans + 1
                                    if ans > L:
                                        return ans
    return ans

def print_L(L):
    m = int(log2(L))
    mc_L = [[]]
    for i in range(m):
        mc_L += [[]]
    if L >= 2:
        # 1 bit
        for i0 in range(L):
            bits_num = [i0]
            ans = leqNum(bits_num, L)
            if ans < L:
                indice = 1 if bits_num[-1] == 0 else ceil(log2(bits_num[-1] + 1))
                mc_L[indice].append(bits_num)
    if L >= 4:
        # 2 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                bits_num = [i0, i1]
                ans = leqNum(bits_num, L)
                if ans < L:
                    indice = 1 if bits_num[-1] == 0 else ceil(log2(bits_num[-1] + 1))
                    mc_L[indice].append(bits_num)
    if L >= 8:
        # 3 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                for i2 in range(i1+1, L):
                    bits_num = [i0, i1, i2]
                    ans = leqNum(bits_num, L)
                    if ans < L:
                        indice = 1 if bits_num[-1] == 0 else ceil(log2(bits_num[-1] + 1))
                        mc_L[indice].append(bits_num)
    if L >= 16:
        # 4 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                for i2 in range(i1+1, L):
                    for i3 in range(i2+1, L):
                        bits_num = [i0, i1, i2, i3]
                        ans = leqNum(bits_num, L)
                        if ans < L:
                            indice = 1 if bits_num[-1] == 0 else ceil(log2(bits_num[-1] + 1))
                            mc_L[indice].append(bits_num)
    if L >= 32:
        # 5 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                for i2 in range(i1+1, L):
                    for i3 in range(i2+1, L):
                        for i4 in range(i3+1, L):
                            bits_num = [i0, i1, i2, i3, i4]
                            ans = leqNum(bits_num, L)
                            if ans < L:
                                indice = 1 if bits_num[-1] == 0 else ceil(log2(bits_num[-1] + 1))
                                mc_L[indice].append(bits_num)
    if L >= 64:
        # 6 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                for i2 in range(i1+1, L):
                    for i3 in range(i2+1, L):
                        for i4 in range(i3+1, L):
                            for i5 in range(i4+1, L):
                                bits_num = [i0, i1, i2, i3, i4, i5]
                                ans = leqNum(bits_num, L)
                                if ans < L:
                                    indice = 1 if bits_num[-1] == 0 else ceil(log2(bits_num[-1] + 1))
                                    mc_L[indice].append(bits_num)
    print_mc_vec(mc_L, L)
    print_mc_flip(mc_L, L)


def print_mc_vec(mc_L, L):
    for s in mc_L:
        print(s)


def print_mc_flip(mc_L, L):
    for s in mc_L:
        print(s)



L = 8
print_L(L)

