"""
运行第一章算法的函数
"""

import Chapter1

petrophila_ori = """aactctatacctcctttttgtcgaatttgtgtgatttatagagaaaatcttattaactga
aactaaaatggtaggtttggtggtaggttttgtgtacattttgtagtatctgatttttaa
ttacataccgtatattgtattaaattgacgaacaattgcatggaattgaatatatgcaaa
acaaacctaccaccaaactctgtattgaccattttaggacaacttcagggtggtaggttt
ctgaagctctcatcaatagactattttagtctttacaaacaatattaccgttcagattca
agattctacaacgctgttttaatgggcgttgcagaaaacttaccacctaaaatccagtat
ccaagccgatttcagagaaacctaccacttacctaccacttacctaccacccgggtggta
agttgcagacattattaaaaacctcatcagaagcttgttcaaaaatttcaatactcgaaa
cctaccacctgcgtcccctattatttactactactaataatagcagtataattgatctga"""

# petrophila_ori_table = Chapter1.frequent_table(text=petrophila_ori, k=9)
# print(petrophila_ori_table)

# 用大肠杆菌 DNA 序列测试 find_clumps 函数
# f = open("E_coli.txt", "r")
# e_coli = f.readline().strip()
# print(len(e_coli))
# e_coli_clumps = Chapter1.find_clumps(e_coli, 9, 500, 5)
# print(e_coli_clumps)

# 测试 approximate_pattern_count 函数
# test_text = "AACAAGCTGATAAACATTTAAAGAG"
# test_pattern = "AAAAA"
# print(Chapter1.approximate_pattern_count(test_text, test_pattern, 1))

# 测试 neighbors 函数
# pattern = "ACGCA"
# test_nbs = Chapter1.neighbors(pattern, 1)
# print(len(test_nbs))
# print(test_nbs)

# 测试 frequent_words_with_mismatches 函数
# test_text = "AACAAGCTGATAAACATTTAAAGAG"
# print(Chapter1.frequent_words_with_mismatches(test_text, 5, 2))

# 测试 frequent_words_with_mismatches_rc 函数
# e_coli = """aatgatgatgacgtcaaaaggatccggataaaacatggtgattgcctcgcataacgcggtatgaaaatggattgaagcccgggccgtggattctactcaactttgtcggcttgagaaagacctgggatcctgggtattaaaaagaagatctatttatttagagatctgttctattgtgatctcttattaggatcgcactgcccTGTGGATAAcaaggatccggcttttaagatcaacaacctggaaagga
# tcattaactgtgaatgatcggtgatcctggaccgtataagctgggatcag
# aatgaggggTTATACACAactcaaaaactgaacaacagttgttcTTTGGA
# TAActaccggttgatccaagcttcctgacagagTTATCCACAgtagatcg
# cacgatctgtatacttatttgagtaaattaacccacgatcccagccattc
# ttctgccggatcttccggaatgtcgtgatcaagaatgttgatcttcagtg""".upper().replace("\n", "")
# print(Chapter1.frequent_words_with_mismatches_rc(e_coli, 9, 1))