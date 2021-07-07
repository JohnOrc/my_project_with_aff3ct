# -*- coding: utf-8 -*-
# @Time    : 2020/12/17 9:35
# @Author  : zyl
# @File    : smallestNumsEven.py
# @Contact : hit_zyl@126.com
# -*- coding: utf-8 -*-
# @Time    : 2020/12/9 15:20
# @Author  : zyl
# @File    : smallestNums.py
# @Contact : hit_zyl@126.com
import numpy as np


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

def print_64(print_type, N):
    L = 64
    # m = np.int(np.log2(L)) = 4
    count = 1
    if print_type == 0:
        print('{}', '\t', 0)
    elif print_type == 1:
        print('pm({})={};'.format(count, 0))
    else:
        print(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    count = 2
    # 2 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            ans = leqNum([i0, i1], L)
            if ans < L:
                if print_type == 0:
                    print('{', i0 + 1, i1 + 1, '}', '\t', ans)
                elif print_type == 1:
                    print('pm({})=pm_({})+pm_({});'.format(count, i0+1, i1+1))
                else:
                    p_str = np.ones([L], np.int8)
                    p_str[i0], p_str[i1] = -1, -1
                    print(' '.join(str(p_str[i]) for i in range(L)))
                count = count + 1
    # 4 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            for i2 in range(i1+1, L):
                for i3 in range(i2+1, L):
                    ans = leqNum([i0, i1, i2, i3], L)
                    if ans < L:
                        
                        print_code(print_type, bits_num)
                        
                        if print_type == 0:
                            print('{', i0 + 1, i1 + 1, i2 + 1, i3 + 1, '}', '\t', ans)
                        elif print_type == 1:
                            print('pm({})=pm_({})+pm_({})+pm_({})+pm_({});'.format(count, i0+1, i1+1, i2+1, i3+1))
                        else:
                            p_str = np.ones([L], np.int8)
                            p_str[i0], p_str[i1], p_str[i2], p_str[i3] = -1, -1, -1, -1
                            print(' '.join(str(p_str[i]) for i in range(L)))
                        count = count + 1
    # 6 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            for i2 in range(i1+1, L):
                for i3 in range(i2+1, L):
                    for i4 in range(i3+1, L):
                        for i5 in range(i4+1, L):
                            ans = leqNum([i0, i1, i2, i3, i4, i5], L)
                            if ans < L:
                                if print_type == 0:
                                    print('{', i0 + 1, i1 + 1, i2 + 1, i3 + 1, i4 + 1, i5 + 1, '}', '\t', ans)
                                elif print_type == 1:
                                    print('pm({})=pm_({})+pm_({})+pm_({})+pm_({})+pm_({})+pm_({});'.format(count, i0+1, i1+1, i2+1, i3+1, i4+1, i5+1))
                                else:
                                    p_str = np.ones([L], np.int8)
                                    p_str[i0], p_str[i1], p_str[i2], p_str[i3], p_str[i4], p_str[i5] = -1, -1, -1, -1, -1, -1
                                    print(' '.join(str(p_str[i]) for i in range(L)))
                                count = count + 1


def print_32(print_type, N):
    L = 32
    # m = np.int(np.log2(L)) = 4
    count = 1
    if print_type == 0:
        print('{}', '\t', 0)
    elif print_type == 1:
        print('pm({})={};'.format(count, 0))
    else:
        print(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    count = 2
    # 2 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            ans = leqNum([i0, i1], L)
            if ans < L:
                if print_type == 0:
                    print('{', i0 + 1, i1 + 1, '}', '\t', ans)
                elif print_type == 1:
                    print('pm({})=pm_({})+pm_({});'.format(count, i0+1, i1+1))
                else:
                    p_str = np.ones([L], np.int8)
                    p_str[i0], p_str[i1] = -1, -1
                    print(' '.join(str(p_str[i]) for i in range(L)))
                count = count + 1
    # 4 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            for i2 in range(i1+1, L):
                for i3 in range(i2+1, L):
                    ans = leqNum([i0, i1, i2, i3], L)
                    if ans < L:
                        if print_type == 0:
                            print('{', i0 + 1, i1 + 1, i2 + 1, i3 + 1, '}', '\t', ans)
                        elif print_type == 1:
                            print('pm({})=pm_({})+pm_({})+pm_({})+pm_({});'.format(count, i0+1, i1+1, i2+1, i3+1))
                        else:
                            p_str = np.ones([L], np.int8)
                            p_str[i0], p_str[i1], p_str[i2], p_str[i3] = -1, -1, -1, -1
                            print(' '.join(str(p_str[i]) for i in range(L)))
                        count = count + 1
    # 6 bits
    if print_type == 0:
        print('{', 1, 2, 3, 4, 5, 6, '}', '\t', L-1)
        print(count)
    elif print_type == 1:
        print('pm({})=pm_({})+pm_({})+pm_({})+pm_({})+pm_({})+pm_({});'.format(count, 1, 2, 3, 4, 5, 6))
    else:
        print(-1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)


def print_16(print_type, N):
    L = 16
    # m = np.int(np.log2(L)) = 4
    count = 1
    if print_type == 0:
        print('{}', '\t', 0)
    elif print_type == 1:
        print('pm({})={};'.format(count, 0))
    else:
        print(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    count = 2
    # 2 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            ans = leqNum([i0, i1], L)
            if ans < L:
                if print_type == 0:
                    print('{', i0 + 1, i1 + 1, '}', '\t', ans)
                elif print_type == 1:
                    print('pm({})=pm_({})+pm_({});'.format(count, i0+1, i1+1))
                else:
                    p_str = np.ones([L], np.int8)
                    p_str[i0], p_str[i1] = -1, -1
                    print(' '.join(str(p_str[i]) for i in range(L)))
                count = count + 1
    # 3 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            for i2 in range(i1+1, L):
                ans = leqNum([i0, i1], L)
                if ans < L:
                    if print_type == 0:
                        print('{', i0 + 1, i1 + 1, i2 + 1, '}', '\t', ans)
                    elif print_type == 1:
                        print('pm({})=pm_({})+pm_({});'.format(count, i0+1, i1+1, i2+1))
                    else:
                        p_str = np.ones([L], np.int8)
                        p_str[i0], p_str[i1], p_str[i2] = -1, -1, -1
                        print(' '.join(str(p_str[i]) for i in range(L)))
                    count = count + 1
    # 4 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            for i2 in range(i1+1, L):
                for i3 in range(i2+1, L):
                    ans = leqNum([i0, i1, i2, i3], L)
                    if ans < L:
                        if print_type == 0:
                            print('{', i0 + 1, i1 + 1, i2 + 1, i3 + 1, '}', '\t', ans)
                        elif print_type == 1:
                            print('pm({})=pm_({})+pm_({})+pm_({})+pm_({});'.format(count, i0+1, i1+1, i2+1, i3+1))
                        else:
                            p_str = np.ones([L], np.int8)
                            p_str[i0], p_str[i1], p_str[i2], p_str[i3] = -1, -1, -1, -1
                            print(' '.join(str(p_str[i]) for i in range(L)))
                        count = count + 1


def print_8(print_type, N):
    # 不考虑0 0-L-1~1-L
    L = 8

    # 0 bit
    bits_num = [ ]
    print_code(print_type, bits_num)
    # 1 bit
    for i0 in range(L):
        ans = leqNum([i0], L)
        if ans < L:
            bits_num = [i0]
            print_code(print_type, bits_num)
    # 2 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            ans = leqNum([i0, i1], L)
            if ans < L:
                bits_num = [i0, i1]
                print_code(print_type, bits_num)
    # 3 bits
    for i0 in range(L):
        for i1 in range(i0+1, L):
            for i2 in range(i1+1, L):
                ans = leqNum([i0, i1, i2], L)
                if ans < L:
                    bits_num = [i0, i1, i2]
                    print_code(print_type, bits_num)

def print_code(print_type, bits_nums):
    l = len(bits_nums)
    print(bits_nums)




print_type = 1
N = 256
# print_4(print_type, N)
# print_8(print_type, N)
# print_16(print_type, N)
# print_32(print_type, N)
print_8(print_type, N)

