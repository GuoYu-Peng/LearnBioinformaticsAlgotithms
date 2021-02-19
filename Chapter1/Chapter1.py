"""
第一章算法
"""


def pattern_count(text: str, pattern: str) -> int:
    """
    给定序列 text 和 k-mer pattern 返回 k-mer 出现次数。
    注意 text 的 index 是 0 base
    :param text: 序列
    :param pattern: k-mer 序列
    :return: k-mer 计数
    """
    kmer_count = 0
    kmer_len = len(pattern)
    k_last = len(text) - kmer_len
    # 因为 range 函数不包含右边界需要加 1
    for i in range(0, k_last + 1):
        text_i = text[i: i + kmer_len]
        if text_i == pattern:
            kmer_count += 1
        else:
            continue
    return kmer_count


# 测试 pattern_count 函数
# text1 = "AACTAAACCTGTGA"
# pattern1 = "AA"
# print(pattern_count(text1, pattern1))


def frequent_words(text: str, k: int) -> dict:
    """
    给定序列 text 和 k-mer 长度 k 给出频率最高的 k-mer 序列
    :param text: 序列
    :param k: k-mer 长度
    :return: 频率最高 k-mer
    """
    kmer_last = len(text) - k
    count_list = list()
    freq_pattern = dict()
    for i in range(0, kmer_last + 1):
        pattern_i = text[i: (i + k)]
        count_i = pattern_count(text, pattern_i)
        count_list.insert(i, count_i)
    max_count = max(count_list)
    max_count_index = list()
    for i, v in enumerate(count_list):
        if v == max_count:
            max_count_index.append(i)
        else:
            continue
    # 利用字典的 key 唯一性实现去重复
    for i in range(0, len(max_count_index)):
        count_index = max_count_index[i]
        pattern_i = text[count_index: (count_index + k)]
        count_i = pattern_count(text, pattern_i)
        freq_pattern[pattern_i] = count_i
    print(freq_pattern)
    return freq_pattern


# 测试 frequent_words 函数
# text2 = "gaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtga"
# frequent_words(text2, 3)


def frequent_table(text: str, k: int) -> dict:
    """
    给定序列 text 和 k-mer 长度 k 给出所有 k-mer 的频数表储存在字典。
    和 frequent_words 函数相比只要遍历一次 text 就能完成，显著降低了计算量。
    :param text: 序列
    :param k: k-mer 长度 k
    :return: k-mer 频数字典
    """
    kmer_freq = dict()
    kmer_last = len(text) - k
    for i in range(0, kmer_last + 1):
        pattern_i = text[i: i + k]
        if pattern_i in kmer_freq.keys():
            kmer_freq[pattern_i] += 1
        else:
            kmer_freq[pattern_i] = 1
    return kmer_freq


# 测试 frequent_table 函数
# text3 = "gaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtga"
# frequent_table(text3, 3)


def better_frequent_words(text: str, k: int) -> dict:
    """
    利用 frequent_table 函数实现计算量更少的 frequent_words 函数改进。
    :param text: 序列
    :param k: k-mer 长度 k
    :return: 最高频 k-mer, 存储在字典
    """
    freq_pattern = dict()
    kmer_freq = frequent_table(text, k)
    max_count = max(kmer_freq.values())
    for k, v in kmer_freq.items():
        if v == max_count:
            freq_pattern[k] = v
        else:
            continue
    print(freq_pattern)
    return freq_pattern


# 测试 better_frequent_words 函数
# better_frequent_words(text2, 3)

# 用 Vibrio_cholerae 的基因组进行测试
# txt = open("Vibrio_cholerae.txt", "r")
# vc_ref = txt.readline().strip()
# better_frequent_words(vc_ref, 9)

# 用 Vibrio cholerae 复制起点区域测试
# vc_ori = """atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaac
# ctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgacca
# cggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgactt
# gtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggatt
# acgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttagga
# tagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaat
# tgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaag
# atcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtt
# tccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"""
# better_frequent_words(vc_ori, 9)


def reverse_complement(text: str) -> str:
    """
    求输入 DNA 序列的反向互补序列
    :param text: DNA 序列
    :return: 反向互补序列
    """
    base_map = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rc_dna = []
    if len(text) == 0:
        print("输入 DNA 序列至少要求 1 碱基")
        rc_dna2 = text
    else:
        for i in range(0, len(text)):
            base = text[i]
            rc_base = base_map[base]
            rc_dna.insert(0, rc_base)
        rc_dna2 = "".join(rc_dna)
    return rc_dna2


def pattern_occurrence(genome: str, pattern: str) -> list:
    """
    给定基因组序列 genome 和子序列 pattern, 寻找出所有子序列位置
    :param genome: 基因组序列
    :param pattern: 子序列
    :return: 子序列出现位置，列表
    """
    last_kmer = len(genome) - len(pattern)
    pattern_position = []
    for i in range(0, last_kmer + 1):
        kmer = genome[i: i + len(pattern)]
        if kmer == pattern:
            pattern_position.append(i)
        else:
            continue
    return pattern_position


# 用 Vibrio_cholerae 的基因组进行测试
# txt = open("Vibrio_cholerae.txt", "r")
# vc_ref = txt.readline().strip()
# print(pattern_occurrence(vc_ref, "ATGATCAAG"))


def find_clumps(text: str, k: int, window_size: int, t: int) -> list:
    """
    从基因组 text 给定窗口长度，kmer 长度，和最小频次 t. 寻找在窗口长度内出现频次超过 t 的 kmer 序列。
    :param text: 基因组序列
    :param k: kmer 长度
    :param window_size: 滑窗长度
    :param t: 最低频次
    :return: kmer 和对应 clumps
    """
    clumps = list()
    window_number = len(text) - window_size + 1
    for w in range(0, window_number):
        window_text = text[w: w + window_size]
        window_table = frequent_table(window_text, k)
        for kmer, freq in window_table.items():
            if freq >= t:
                clumps.append(kmer)
            else:
                continue
    kmer_clump = list(set(clumps))
    return kmer_clump


def base_score(base: str) -> int:
    if base == "G":
        return 1
    elif base == "C":
        return -1
    else:
        return 0


def minimum_skew(sequence: str) -> list:
    """
    Skew 定义为从 0-i 长度的 DNA 序列 G, C 碱基数目差，这个函数找到给定序列 Skew 最小的位置。
    :param sequence: DNA 序列
    :return: 最小 Skew 位置
    """
    seq_len = len(sequence)
    skew = [0 for x in range(0, seq_len + 1)]
    base_1 = sequence[0]
    skew[1] = base_score(base_1)
    for i in range(2, seq_len):
        base_i = sequence[i]
        score_i = base_score(base_i)
        skew[i] = skew[i - 1] + score_i
    min_skew = min(skew)
    min_position = list()
    for i, v in enumerate(skew):
        if v == min_skew:
            min_position.append(i)
        else:
            continue
    return min_position


def hamming_distance(seq_1: str, seq_2: str) -> int:
    """
    计算 2 条等长序列的 hamming distance.
    :param seq_1: DNA 序列 1
    :param seq_2: DNA 序列 2
    :return: Hamming distance
    """
    hd = 0
    seq_length1 = len(seq_1)
    seq_length2 = len(seq_2)
    if seq_length1 != seq_length2:
        print("输入的序列长度不相等")
        exit()
    else:
        for i in range(0, seq_length1):
            p = seq_1[i]
            q = seq_2[i]
            if p != q:
                hd += 1
            else:
                continue
    return hd


def approximate_pattern_count(text: str, pattern: str, d: int) -> int:
    """
    给定序列 text 和 pattern 求在 hamming distance 小于 d 的匹配次数。
    :param text: 长序列
    :param pattern: 目标 pattern
    :param d: 允许的最大 Hamming distance
    :return: 匹配次数
    """
    text_length = len(text)
    k = len(pattern)
    last_kmer_pos = text_length - k
    pattern_counts = 0
    for i in range(0, last_kmer_pos + 1):
        kmer_i = text[i: i + k]
        hd_i = hamming_distance(pattern, kmer_i)
        if hd_i <= d:
            pattern_counts += 1
    return pattern_counts


def one_mutate(sequence: str) -> list:
    """
    给一个序列，输出有一个点突变的序列集合。
    :param sequence: 输入序列
    :return: 有一个点突变的序列集合
    """
    bases = ("A", "T", "C", "G")
    seq_len = len(sequence)
    seq_list = list(sequence)
    mutate_seqs = list()
    for i in range(0, seq_len):
        seq_list2 = seq_list.copy()
        for base in bases:
            if seq_list2[i] == base:
                continue
            else:
                seq_list2[i] = base
                mutate_seq = "".join(seq_list2)
                mutate_seqs.append(mutate_seq)
    return mutate_seqs


def neighbors(pattern: str, d: int) -> list:
    """
    输入序列 pattern 生成所有 hamming distance 在 d 以内的序列
    :param pattern: 输入序列
    :param d: 最大 Hamming distance
    :return: 符合要求的序列集合
    """
    all_seqs = [0 for x in range(0, d)]
    seqs_0 = one_mutate(pattern)
    all_seqs[0] = seqs_0
    for i in range(1, d):
        seqs_i = list()
        for each in all_seqs[i - 1]:
            seqs_i.extend(one_mutate(each))
        all_seqs[i] = seqs_i
    all_seqs2 = list()
    for seq_list in all_seqs:
        all_seqs2.extend(seq_list)
    all_seqs3 = list(set(all_seqs2))
    return all_seqs3


def frequent_words_with_mismatches(text: str, k: int, d: int) -> dict:
    """
    给定序列 text 和 k-mer 长度，分析在错配 d 以内的最高频数 k-mer.
    :param text: 输入序列 text
    :param k: k-mer 长度
    :param d: 最大 Hamming distance
    :return: 频数最高的 k-mer 及频数
    """
    k_mers = list()
    last_position = len(text) - k
    for p in range(0, last_position):
        pattern_p = text[p: p + k]
        k_mers_p = neighbors(pattern_p, d)
        k_mers.extend(k_mers_p)
    k_mers2 = set(k_mers)
    k_mers_counts = dict()
    for k_mer in k_mers2:
        k_mer_count = approximate_pattern_count(text, k_mer, d)
        k_mers_counts[k_mer] = k_mer_count
    max_count = 0
    for k, v in k_mers_counts.items():
        if v > max_count:
            max_count = v
        else:
            continue
    k_mers_counts2 = dict()
    for k, v in k_mers_counts.items():
        if v == max_count:
            k_mers_counts2[k] = v
        else:
            continue
    return k_mers_counts2


def frequent_words_with_mismatches2(text: str, k: int, d: int) -> dict:
    """
    给定序列 text 和 k-mer 长度，分析在错配 d 以内的 k-mer 频数.
    和上面函数区别是不求出最大频数，给出全部结果。
    :param text: 输入序列 text
    :param k: k-mer 长度
    :param d: 最大 Hamming distance
    :return: 频数最高的 k-mer 及频数
    """
    k_mers = list()
    last_position = len(text) - k
    for p in range(0, last_position):
        pattern_p = text[p: p + k]
        k_mers_p = neighbors(pattern_p, d)
        k_mers.extend(k_mers_p)
    k_mers2 = set(k_mers)
    k_mers_counts = dict()
    for k_mer in k_mers2:
        k_mer_count = approximate_pattern_count(text, k_mer, d)
        k_mers_counts[k_mer] = k_mer_count
    return k_mers_counts


def frequent_words_with_mismatches_rc(text: str, k: int, d: int) -> dict:
    """
    给定序列 text 和 k-mer 长度，分析在错配 d 以内的最高频数 k-mer.
    和上面函数相比这个函数考虑了反向互补序列。
    :param text: 输入序列 text
    :param k: k-mer 长度
    :param d: 最大 Hamming distance
    :return: 频数最高的 k-mer 及频数
    """
    rc_text = reverse_complement(text)
    k_mers_count1 = frequent_words_with_mismatches2(text, k, d)
    k_mers_count2 = frequent_words_with_mismatches2(rc_text, k, d)
    k_mers_count3 = dict()
    for k, v in k_mers_count1.items():
        rc_k = reverse_complement(k)
        v2 = k_mers_count2.get(rc_k)
        if v2:
            v3 = v + v2
        else:
            v3 = v
        k_mers_count3[k] = v3
    k_mers1 = k_mers_count1.keys()
    for k, v in k_mers_count2.items():
        rc_k = reverse_complement(k)
        if rc_k in k_mers1:
            continue
        else:
            k_mers_count3[rc_k] = v
    max_count = 0
    for k, v in k_mers_count3.items():
        if v > max_count:
            max_count = v
        else:
            continue
    k_mers_counts4 = dict()
    for k, v in k_mers_count3.items():
        if v == max_count:
            k_mers_counts4[k] = v
        else:
            continue
    return k_mers_counts4