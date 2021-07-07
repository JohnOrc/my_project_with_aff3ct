# -*- coding: utf-8 -*-
# @Time    : 2020/12/17 9:35
# @Author  : zyl
# @File    : mc_generate.py
# @Contact : hit_zyl@126.com
from math import log2, ceil
from argparse import ArgumentParser


def leqNumEven(l_list, L):
    ans = 0
    # 2 bits
    if len(l_list) >= 2:
        for i0 in range(L):
            for i1 in range(i0 + 1, L):
                if i0 <= l_list[-2] and i1 <= l_list[-1]:
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

def print_L(L, f):
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
        print('L = 2')
    if L >= 4:
        # 2 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                bits_num = [i0, i1]
                ans = leqNum(bits_num, L)
                if ans < L:
                    indice = ceil(log2(bits_num[-1] + 1))
                    mc_L[indice].append(bits_num)
        print('L = 4')                
    if L >= 8:
        # 3 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                for i2 in range(i1+1, L):
                    bits_num = [i0, i1, i2]
                    ans = leqNum(bits_num, L)
                    if ans < L:
                        indice = ceil(log2(bits_num[-1] + 1))
                        mc_L[indice].append(bits_num)
        print('L = 8')
    if L >= 16:
        # 4 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                for i2 in range(i1+1, L):
                    for i3 in range(i2+1, L):
                        bits_num = [i0, i1, i2, i3]
                        ans = leqNum(bits_num, L)
                        if ans < L:
                            indice = ceil(log2(bits_num[-1] + 1))
                            mc_L[indice].append(bits_num)
        print('L = 16')
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
                                indice = ceil(log2(bits_num[-1] + 1))
                                mc_L[indice].append(bits_num)
        print('L = 32')
    if L == 64:
        mc_L[m].append([0,1,2,3,4,5])
        print('L = 64')

    print_mc_vec(mc_L, L, f)

    f.write('\n\n\n')

    print_mc_flip(mc_L, L, f)

def print_L_even(L, f):
    m = int(log2(L))
    mc_L = [[]]
    for i in range(m):
        mc_L += [[]]

    if L == 2:
        mc_L[0].append([0,1])
        print('L = 2')

    if L >= 4:
        # 2 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                bits_num = [i0, i1]
                ans = leqNumEven(bits_num, L)
                if ans < L:
                    indice = ceil(log2(bits_num[-1] + 1))
                    mc_L[indice].append(bits_num)
        print('L = 4')    

    if L == 8:
        # 4 bits
        bits_num = [0,1,2,3]
        indice = ceil(log2(bits_num[-1] + 1))
        mc_L[indice].append(bits_num)
        print('L = 8')

    if L >= 16:
        # 4 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                for i2 in range(i1+1, L):
                    for i3 in range(i2+1, L):
                        bits_num = [i0, i1, i2, i3]
                        ans = leqNumEven(bits_num, L)
                        if ans < L:
                            indice = ceil(log2(bits_num[-1] + 1))
                            mc_L[indice].append(bits_num)
        print('L = 16')

    if L == 32:
        # 6 bits        
        bits_num = [0,1,2,3,4,5]
        indice = ceil(log2(bits_num[-1] + 1))
        mc_L[indice].append(bits_num)
        print('L = 32')

    flag = False
    if L >= 64:
        # 6 bits
        for i0 in range(L):
            for i1 in range(i0+1, L):
                for i2 in range(i1+1, L):
                    for i3 in range(i2+1, L):
                        for i4 in range(i3+1, L):
                            for i5 in range(i4+1, L):
                                bits_num = [i0, i1, i2, i3, i4, i5]
                                ans = leqNumEven(bits_num, L)
                                if ans < L:
                                    indice = ceil(log2(bits_num[-1] + 1))
                                    mc_L[indice].append(bits_num)
                                    if bits_num[0] == 1:
                                        flag = True
                                        break
                                else:
                                    break
                            if flag:
                                break
                        if flag:
                            break
                    if flag:
                        break
                if flag:
                    break
            if flag:
                break
        print('L = 64')
        

    print_mc_vec_even(mc_L, L, f)

    f.write('\n\n\n')

    print_mc_flip_even(mc_L, L, f)

def find_num(vec_list, vec):
    for i, v in enumerate(vec_list):
        if v == vec:
            return i
    return len(vec_list)


def print_mc_vec_even(mc_L, L, f):
    count = 0
    vec_list = ['e']


    for n_stage in mc_L:
        for bits in n_stage:
            count = count + 1
            vec_list.append('e')
            for b in bits:
                vec_list[count] = vec_list[count] + '+a{}'.format(b)


    f.write('metrics_vec[1][c_num * path + 0] = metrics [path]; // empty\n')
    count = 0
    tab = ''
    for n, n_stage in enumerate(mc_L):
        if n >= 2:
            f.write('if (n_elmts >= {})\n'.format(1 << (n)))
            f.write('{\n')
            tab = tab + '\t'
            
        for bits in n_stage:
            s = 'e'
            for b in bits[:-2]:
                s = s + '+a{}'.format(b)
            num = find_num(vec_list, s)
            count = count + 1
            f.write(tab + 'metrics_vec[1][c_num * path + {}] = sat_m<R>(metrics_vec[1][c_num * path + {}] + pen[{}] + pen[{}]);'.format(count, num, bits[-2], bits[-1]) + ' // ' + s + '+a{}+a{}\n'.format(bits[-2], bits[-1]))
            



def print_mc_flip_even(mc_L, L, f):
    f.write('switch( dup )\n')
    f.write('{\n')
    f.write('case 0:\n')
    f.write('\t// nothing to do\n')
    f.write('\tbreak;\n')

    count = 1
    for n_stage in mc_L:
        for bits in n_stage:
            s = '// e '
            for b in bits:
                s = s + '+ a{} '.format(b)
            f.write('case {}:'.format(count) + s + '\n')
            count = count + 1
            for b in bits:
                f.write('\ts[new_path][off_s + bit_flips_r1[bits_num * old_path +{}]] = s[old_path][off_s + bit_flips_r1[bits_num * old_path +{}]] ? 0 : b;\n'.format(b, b))
                f.write('\tbreak;\n')
    f.write('default:\n')
    f.write('\tthrow tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on rate 1 node.");\n')
    f.write('\tbreak;\n')
    f.write('}')


def print_mc_vec(mc_L, L, f):
    count = 0
    vec_list = ['e']
    for n_stage in mc_L:
        for bits in n_stage:
            count = count + 1
            vec_list.append('e')
            for b in bits:
                vec_list[count] = vec_list[count] + '+a{}'.format(b)

    f.write('metrics_vec[1][c_num * path + 0] = metrics [path]; // empty\n')
    count = 0
    tab = ''
    for n, n_stage in enumerate(mc_L):
        if n >= 2:
            f.write('if (n_elmts >= {})\n'.format(1 << (n)))
            f.write('{\n')
            tab = tab + '\t'
            
        for bits in n_stage:
            s = 'e'
            for b in bits[:-1]:
                s = s + '+a{}'.format(b)
            num = find_num(vec_list, s)
            count = count + 1
            f.write(tab + 'metrics_vec[1][c_num * path + {}] = sat_m<R>(metrics_vec[1][c_num * path + {}] + pen[{}]);'.format(count, num, bits[-1]) + ' // ' + s + '+a{}\n'.format(bits[-1]))
            



def print_mc_flip(mc_L, L, f):
    f.write('switch( dup )\n')
    f.write('{\n')
    f.write('case 0:\n')
    f.write('\t// nothing to do\n')
    f.write('\tbreak;\n')

    count = 1
    for n_stage in mc_L:
        for bits in n_stage:
            s = '// e '
            for b in bits:
                s = s + '+ a{} '.format(b)
            f.write('case {}:'.format(count) + s + '\n')
            count = count + 1
            for b in bits:
                f.write('\ts[new_path][off_s + bit_flips_r1[bits_num * old_path +{}]] = s[old_path][off_s + bit_flips_r1[bits_num * old_path +{}]] ? 0 : b;\n'.format(b, b))
                f.write('\tbreak;\n')
    f.write('default:\n')
    f.write('\tthrow tools::runtime_error(__FILE__, __LINE__, __func__, "Flip bits error on rate 1 node.");\n')
    f.write('\tbreak;\n')
    f.write('}')


def main():
    L = 8
    parser = ArgumentParser()
    parser.add_argument('--l', type=int, dest='list_size', metavar='L', default=L)
    params = parser.parse_args()


    f = open('./{}even.txt'.format(params.list_size), 'w')

    # print_L(L=params.list_size, f=f)
    print_L_even(L=params.list_size, f=f)

    f.close()

if __name__ == '__main__':
    main()