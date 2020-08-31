#!/usr/bin/python3

import spot
import enum
import sys
import argparse
from scipy.optimize import linear_sum_assignment

class MarkType(enum.Enum):
    """Represents whether acceptance set T appears in the formula as Fin(T) or Inf(T).
    
    Arguments:
        enum {int} 
    """    
    Inf = 1
    Fin = 2


class ACCMark: 
    """Represents an acceptance set.
    """        
    def __init__(self, mtype, num):
        """ACCMark constructor.
        
        Arguments:
            mtype {MarkType} -- denotes whether a set appears in Inf or Fin term in the formula
            num {int} -- number of the acceptance set
        """        
        self.type = mtype
        self.num = num

    def __str__(self):
        """Return the information about ACCMark as string.
        
        Returns:
            string -- ACCMarks information as string
        """        
        return ''.join(["Inf" if self.type == MarkType.Inf else "Fin", '(', str(self.num), ')'])

    def __eq__(self, other): 
        """Compare two ACCMarks, return true if they are the same, false otherwise.
        
        Arguments:
            other {ACCMark} -- the other ACCMark
        
        Returns:
            bool -- true if marks are the same, false otherwise
        """        
        if (self.type == other.type) and (self.num == other.num): 
            return True
        else: 
            return False


class PACC: 
    def __init__(self, acc):
        self.formula = parse_acc(acc)
        self.sat = None

    def __str__(self):
        f = []
        for dis in self.formula:
            f.append('(')
            for con in dis:
                f.append(str(con))
                if con is not dis[-1]:
                    f.append(" & ")
            f.append(')')
            if dis is not self.formula[-1]:
                f.append(" | ")
        return ''.join(f)

    def __getitem__(self, index):
        return self.formula[index]

    def __len__(self):
        return len(self.formula)

    def acc_len(self):
        l = 0
        for dis in self.formula:
            for con in dis:
                l += 1
        return l

    def set_sat(self, val):
        self.sat = val

    def count_unique_m(self):                
        inf = []
        fin = []
        for dis in self.formula:
            for con in dis:
                if con.type == MarkType.Inf and con.num not in inf:
                    inf.append(con.num)
                if con.type == MarkType.Fin and con.num not in fin:
                    fin.append(con.num)
        return len(inf), len(fin)

    def count_total_unique_m(self):                
        um = []
        for dis in self.formula:
            for con in dis:
                if con not in um:
                    um.append(con)              
        return len(um)

    def max(self):
        m = self.formula[0][0].num
        for dis in self.formula:
            for con in dis:
                if con.num > m:
                    m = con.num
        return m
                
    def int_format(self):
        f = []
        for dis in self.formula:
            d = []
            for con in dis:
                d.append(con.num)
            f.append(d)
        return f

    def resolve_redundancy(self):
        rem_d = []
        int_f = self.int_format()
        for i in range(len(self.formula)):
            for j in range(len(self.formula)):
                if i != j and set(int_f[i]).issubset(set(int_f[j])) and not set(int_f[i]) == set(int_f[j]) and j not in rem_d:
                    rem_d.append(j)

        res_f = []
        for i in range(len(self.formula)):
            if self.formula[i] not in res_f and i not in rem_d:
                res_f.append(self.formula[i])
        self.formula = res_f

    def clean_up(self, aut, scc):
        clean_f = []
        marks = scc_current_marks(aut, scc)
        for dis in self.formula:
            clean_dis = []
            for con in dis:
                if con.num in marks:
                    clean_dis.append(con)
            if clean_dis:
                clean_f.append(clean_dis)
        self.formula = clean_f
        self.resolve_redundancy()

    def get_mtype(self, m):
        for dis in self.formula:
            for con in dis:
                if m == con.num:
                    return con.type

    def find_m_dis(self, m):
        occurrences = []
        i = 0
        for dis in self.formula:
            for con in dis:
                if con.num == m:
                    occurrences.append(i)   
            i += 1   
        return occurrences

    def rem_from_dis(self, index, m):
        new_dis = []
        for con in self.formula[index]:
            if con.num != m:
                new_dis.append(con)
        self.formula[index] = new_dis

    def rem_dis(self, index):
        new_f = []
        i = 0
        for dis in self.formula:
            if i != index:
                new_f.append(dis)
            i += 1
        self.formula = new_f    

    
### PARSE ACC ###

def parse_acc(acc):
    """Parses an acc in DNF and returns the formula represented by list of
    lists of ACCMarks. The inner lists represent disjuncts of the formula. The
    inner lists contain ACCMark (see ACCMark class documentation) objects
    representing atomic conditions (such as Inf(1)).
    
    Arguments:
        acc {spot::acc_cond::acc_code} -- spot representation of acceptance condition formula
    
    Returns:
        [[ACCMark]] -- list of lists of ACCMarks
    """    
    if acc is "": 
        return []
    formula = []
    for dis in str(acc).split('|'):
        new_dis = []
        for con in dis.split('&'):
            mtype = None
            if con.find("Fin") != -1:
                mtype = MarkType.Fin
            else:
                mtype = MarkType.Inf
            num = []
            for c in con:
                if c.isdigit():  
                    num.append(c)
            new_dis.append(ACCMark(mtype, int(''.join(num))))               
        formula.append(new_dis)            
    return formula


### SIMPLIFY AUXILIARY ###

def scc_clean_up_edges(aut, acc, scc):
    """Removes useless marks from the edges of the given SCC. Useless marks denote transitions of sets that do not appear in the Acc.
    
    Arguments:
        aut {spot::twa} -- input automaton
        acc {PACC} -- acceptance condition formula
        scc {spot::scc_info_node} -- SCC
    """    
    current_m = []
    for dis in acc.formula:
        for con in dis:
            if con.num not in current_m:
                current_m.append(con.num)

    for s in scc.states():
        for e in aut.out(s):
            for m in e.acc.sets():
                if int(m) not in current_m:
                    e.acc.clear(m)




def scc_current_marks(aut, scc): 
    """Return a list of marks currently present on edges in given scc.
    
    Arguments:
        aut {spot::twa} -- input automaton
        scc {spot::scc_info_node} -- SCC
    
    Returns:
        [int] -- list of marks on edges of given scc
    """    
    marks = []
    for s in scc.states():
        for e in aut.out(s): 
            for m in e.acc.sets():
                if m not in marks:
                    marks.append(int(m)) 
    return marks


def scc_everywhere(aut, scc):
    """Return a list of acceptance sets that include all transitions in the SCC.
    
    Arguments:
        aut {spot::twa} -- input automaton
        scc {spot::scc_info_node} -- SCC
    
    Returns:
        [int] -- list of marks (acceptance sets) that include all transtitions of SCC
    """    
    m_all_edges = set(scc_current_marks(aut, scc))
    for s in scc.states():
        for e in aut.out(s):
            if e.dst in scc.states():
                for m in list(m_all_edges):
                    if not m in e.acc.sets():
                        m_all_edges.remove(m) 
    return list(m_all_edges)


def scc_compl_sets(aut, scc): 
    """Returns a list of tuples. Each tuple contains numbers of acceptance sets m1, m2, 
    such that the set of transitions that hold m1 is complementary with the set of transitions
    that hold m2.
    
    Returns:
        [(int, int))] -- array of tuples of complementary marks in given scc 
    """    
    compl = []
    marks = scc_current_marks(aut, scc)
    for m1 in marks:
        for m2 in marks:
            if m1 is not m2:            
                are_compl = True
                for s in scc.states():
                    for e in aut.out(s):
                        if e.dst in scc.states() and ((m1 in e.acc.sets() and m2 in e.acc.sets()) or (m1 not in e.acc.sets() and m2 not in e.acc.sets())):
                            are_compl = False
                if are_compl and (m1, m2) not in compl and (m2, m1) not in compl:
                    compl.append((m1, m2))
    return compl


def scc_subsets(aut, scc):
    """[summary]
    
    Arguments:
        aut {spot::twa} -- input automaton
        scc {spot::scc_info_node} -- SCC
    
    Returns:
        [(int, int)] -- list of tuples that denote complementary sets
    """    
    marks = scc_current_marks(aut, scc)
    subsets = []
    for m1 in marks:
        for m2 in marks:
            if m1 is not m2:
                is_sub = True
                for s in scc.states():
                    for e in aut.out(s):
                        if e.dst in scc.states() and m2 in e.acc.sets() and not m1 in e.acc.sets():
                            is_sub = False
                if is_sub:
                    subsets.append((m1, m2))
    return subsets 
                     
                        
def simpl_inf_con(aut, acc, scc, subsets):    
    """If there are acceptance sets X, Y, such that X is a subset of Y, and X, Y always appear 
    in the same disjuncts, then remove X.
    
    Arguments:
        aut {spot::twa} -- input automaton
        acc {PACC} -- acceptance condition formula
        scc {spot::scc_info_node} -- SCC
        subsets {[(int, int)]} -- list if tuples that denote which set includes which
    """     
    for sub in subsets:        
        if acc.get_mtype(sub[0]) == MarkType.Inf and acc.get_mtype(sub[1]) == MarkType.Inf:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            if (all(i in sub_i for i in super_i)):
                for index in reversed(super_i):
                    acc.rem_from_dis(index, sub[0])                
                scc_clean_up_edges(aut, acc, scc)
                simpl_inf_con(aut, acc, scc, scc_subsets(aut, scc))
                return                     


def simpl_fin_same_dis(aut, acc, scc):
    fins = []
    for dis in acc.formula:
        for con in dis:
            if con.type == MarkType.Fin:
                fins.append(con)
                
    for fin1 in fins:
        for fin2 in fins:
            if fin1 is not fin2:
                merge = True
                for dis in acc.formula:
                    if (fin1 not in dis and fin2 in dis) or (fin2 not in dis and fin1 in dis):
                        merge = False
                if merge:
                    add_dupl_marks(aut, scc, fin2.num, fin1.num)
                    acc.clean_up(aut, scc)
                    return


def simpl_fin_con_subsets(aut, acc, scc, subsets):
    for sub in subsets:
        if acc.get_mtype(sub[0]) == MarkType.Fin and acc.get_mtype(sub[1]) == MarkType.Fin:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            if (all(i in super_i for i in sub_i)):
                for index in sub_i:
                    acc.rem_from_dis(index, sub[1])
                scc_clean_up_edges(aut, acc, scc)
                acc.clean_up(aut, scc)
                simpl_fin_con_subsets(aut, acc, scc, scc_subsets(aut, scc))
                return
            else:
                for i in reversed(sub_i):
                    if i in super_i:
                        acc.rem_from_dis(i, sub[1])
                        acc.clean_up(aut, scc)                           


def simpl_co_con(aut, acc, scc, compl_sets): 
    for co in compl_sets:
        if acc.get_mtype(co[0]) != acc.get_mtype(co[1]):
            inf, fin = co[0], co[1] 
            if (acc.get_mtype(co[0]) == MarkType.Fin):
                 inf, fin = co[1], co[0]
            inf_i = acc.find_m_dis(inf)
            fin_i = acc.find_m_dis(fin)
            if all(i in fin_i for i in inf_i):                
                for index in reversed(inf_i):
                    acc.rem_from_dis(index, inf)
                
                scc_clean_up_edges(aut, acc, scc)
                simpl_co_con(aut, acc, scc, scc_compl_sets(aut, scc))
                return
            else:
                for i in inf_i:
                    if i in fin_i:
                        acc.rem_from_dis(i, inf)
                acc.clean_up(aut, scc)


#check if inf marks are placed on all edges that do not have fin mark
def check_implies_marks(aut, scc, fin, inf):
    for s in scc.states():
        for e in aut.out(s):
            if e.dst in scc.states() and not inf in e.acc.sets() and not fin in e.acc.sets():
                return False
    return True


#(Fin AND Inf), situations where: 'if Fin true then Inf true'; Inf can be removed
def simpl_fin_implies_inf(aut, acc, scc):
    for i in range(0, len(acc)):
        infs = []
        for con in acc[i]:
            if con.type == MarkType.Inf:
                infs.append(con)
        if not infs:
            return
        for con in acc[i]:
            if con.type == MarkType.Fin:
                for inf in infs:
                    if check_implies_marks(aut, scc, con.num, inf.num):                        
                        acc.rem_from_dis(i, inf.num)
                        scc_clean_up_edges(aut, acc, scc)
                        simpl_fin_implies_inf(aut, acc, scc)
                        return
    
    
def simpl_substitute(aut, acc, scc):
    current = scc_current_marks(aut, scc)
    everywhere = scc_everywhere(aut, scc)
    rem_d = []
    for i in range(len(acc)):
        for con in acc[i]:
            if (con.num in everywhere and con.type == MarkType.Fin) or (con.num not in current and con.type == MarkType.Inf):
                if i not in rem_d:
                    rem_d.append(i)                
            elif con.num in everywhere and con.type == MarkType.Inf:
                acc.rem_from_dis(i, con.num)
            elif con.num not in current and con.type == MarkType.Fin:
                acc.rem_from_dis(i, con.num)
    for index in reversed(rem_d):
        acc.rem_dis(index)
    acc.clean_up(aut, scc)


def simpl_false_subsets(aut, acc, subsets, scc):
    for sub in subsets:
        if acc.get_mtype(sub[0]) == MarkType.Fin and acc.get_mtype(sub[1]) == MarkType.Inf:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            for i in sub_i:
                if i in super_i:
                    acc.rem_dis(i)
                    scc_clean_up_edges(aut, acc, scc)                    
                    simpl_false_subsets(aut, acc, scc_subsets(aut, scc), scc)
                    return


def simpl_false_fin(aut, acc, scc):
    rem_dis = []
    for i in range(len(acc)):
        current_fins = []
        for con in acc[i]:
            if con.type == MarkType.Fin:
                current_fins.append(con.num)
        current_dis = True
        for s in scc.states():
            for e in aut.out(s):
                if e.dst in scc.states() and all(m not in e.acc.sets() for m in current_fins):
                    current_dis = False
        if current_dis:
            rem_dis.append(i)
    for index in reversed(rem_dis):
        acc.rem_dis(index)


### SIMPLIFY ###

def simplify(aut, acc, scc):
    acc_l = acc.count_total_unique_m()

    # remove always false disjuncts
    simpl_substitute(aut, acc, scc)
    simpl_false_subsets(aut, acc, scc_subsets(aut, scc), scc)
    simpl_false_fin(aut, acc, scc)

    # simplify inclusion mark sets
    simpl_inf_con(aut, acc, scc, scc_subsets(aut, scc))
    simpl_fin_con_subsets(aut, acc, scc, scc_subsets(aut, scc))
    simpl_fin_same_dis(aut, acc, scc)
    
    # simplify complementary marks
    simpl_co_con(aut, acc, scc, scc_compl_sets(aut, scc))

    #simplify disjuncts where (Fin = True) => (Inf = True)
    simpl_fin_implies_inf(aut, acc, scc)

    scc_clean_up_edges(aut, acc, scc)
    acc.clean_up(aut, scc)

    if acc_l > acc.count_total_unique_m():
        simplify(aut, acc, scc)   
    

### MERGE AUXILIARY ###

def shift_fst_acc(aut, acc, scc, m):
    next_m = m + 1
    log = {}

    for dis in acc.formula:
        for con in dis:
            if con.num in log:
                con.num = log[con.num]
            else:
                add_dupl_marks(aut, scc, con.num, next_m)
                log[con.num] = next_m
                con.num = next_m
                next_m += 1
    scc_clean_up_edges(aut, acc, scc)


def count_cost(dis1, dis2):
    inf, fin = 0, 0
    for con in dis1:
        if con.type == MarkType.Inf:
            inf -= 1
        else:
            fin -= 1
    for con in dis2:
        if con.type == MarkType.Inf:
            inf += 1
        else:
            fin += 1
    if inf < 0:
        inf = 0
    if fin < 0:
        fin = 0
    return inf + fin


def make_matrix(acc1, acc2):
    m = []
    for i in range(len(acc1)):
        row = []
        for j in range(len(acc1)):
            if j >= len(acc2):
                row.append(0)
            else:
                row.append(count_cost(acc1[i], acc2[j]))
        m.append(row)
    return m


def replace_marks(aut, scc, origin_m, new_m):
    for s in scc.states():
        for e in aut.out(s):
            if e.dst in scc.states() and origin_m in e.acc.sets():
                e.acc.clear(origin_m)
                if not new_m in e.acc.sets():
                    e.acc.set(new_m)


def partition(accs, sccs, low, high): 
    i = (low - 1) 
    pivot = len(accs[high])
  
    for j in range(low , high):   
        if   len(accs[j]) > pivot:           
            i += 1 
            accs[i], accs[j] = accs[j], accs[i] 
            sccs[i], sccs[j] = sccs[j], sccs[i]
  
    accs[i + 1], accs[high] = accs[high], accs[i + 1] 
    sccs[i + 1], sccs[high] = sccs[high], sccs[i + 1]
    return (i + 1) 
  
def acc_quicksort(accs, sccs, low, high): 
    if low < high: 
   
        pi = partition(accs, sccs, low, high) 
  
        acc_quicksort(accs, sccs, low, pi - 1) 
        acc_quicksort(accs, sccs, pi + 1, high) 


def place_con(dis, con, used, dis_remaining):
    if con in dis and not used[dis.index(con)]:
        return dis.index(con)
    remaining = [] #list of numbers of marks that remain unpaired (whatever comes after the current conjunct in the clause)
    for c in dis_remaining:
        if con == c or remaining:
            remaining.append(c.num) 
    for i in range(len(dis)):
        if dis[i].type == con.type and not used[i] and dis[i].num not in remaining:
            return i         
    for i in range(len(dis)):
        if dis[i].type == con.type and not used[i]:            
            return i
    return None
        

def add_dupl_marks(aut, scc, origin_m, new_m):
    for s in scc.states():
        for e in aut.out(s):
            if e.dst in scc.states() and origin_m in e.acc.sets():
                if not new_m in e.acc.sets():
                    e.acc.set(new_m)


def resolve_depend(aut, acc1, scc1, acc2, scc2):
    m = make_matrix(acc1, acc2)
    row_i, col_i = linear_sum_assignment(m)
    log = {} #used to find dependencies
    l = {} #return log
    for i in range(len(row_i)):
        dis1 = acc1[col_i[i]]
        used = [False] * len(dis1)
        dis2 = []
        if row_i[i] < len(acc2):
            dis2 = acc2[row_i[i]]
        for con in dis2:
            index = place_con(dis1, con, used, dis2)
            if index is not None:
                if con.num in log:
                    if log[con.num] == dis1[index].num:
                        used[index] = True
                    else:                       
                        used[index] = True
                        l[dis1[index].num] = acc1.max() + 1
                        log[con.num] = acc1.max() + 1
                        add_dupl_marks(aut, scc1, dis1[index].num, acc1.max() + 1)
                        dis1[index].num = acc1.max() + 1                                                                      
                else:
                    log[con.num] = dis1[index].num                
                    used[index] = True
        if not dis2:
            for con in dis1:                
                if con.num in log.values():
                    l[con.num] = acc1.max() + 1                    
                    add_dupl_marks(aut, scc1, con.num, acc1.max() + 1)                   
                    con.num = acc1.max() + 1
        for j in range(len(used)):
            if not used[j]:
                if dis1[j].num in log.values():
                    l[dis1[j].num] = acc1.max() + 1
                    add_dupl_marks(aut, scc1, dis1[j].num, acc1.max() + 1)
                    dis1[j].num = acc1.max() + 1
    return (scc1, l)


### MERGE ###

def merge(aut, acc1, scc1, acc2, scc2):
    m = make_matrix(acc1, acc2)
    col_ind, row_ind = linear_sum_assignment(m) 
    log = (scc2, {})
    for i in range(len(row_ind)):
        dis1 = acc1[col_ind[i]]
        used = [False] * len(dis1)
        dis2 = []
        if row_ind[i] < len(acc2):
            dis2 = acc2[row_ind[i]]
        for con in dis2:
            index = place_con(dis1, con, used, dis2)
            if index is None:
                new_num = acc1.max() + 1
                log[1][con.num] = new_num
                if con.type == MarkType.Inf:                                
                    for s in scc1.states():
                        for e in aut.out(s):
                            if e.dst in scc1.states():
                                e.acc.set(new_num) 
                add_dupl_marks(aut, scc2, con.num, new_num)
                dis1.append(ACCMark(con.type, new_num))
                used.append(True)
            else:            
                log[1][con.num] = dis1[index].num
                add_dupl_marks(aut, scc2, con.num, dis1[index].num)
                used[index] = True
    return log


def merge_accs(aut, sccs, accs):
    nempty_accs = []
    nempty_sccs = []

    for i in range(len(accs)):
        if accs[i].formula:
            nempty_accs.append(accs[i])
            nempty_sccs.append(sccs[i])

    if not nempty_accs:
        return None, None
    
    acc_quicksort(nempty_accs, nempty_sccs, 0, len(nempty_accs) - 1)

    merged_f = PACC(str(nempty_accs[0]))
    for i in range(len(accs)):
        if accs[i] is nempty_accs[0]:
            accs[i] = merged_f
            break
    
    m = 0
    for acc in nempty_accs:
        if acc.max() > m:
            m = acc.max()
    shift_fst_acc(aut, merged_f, nempty_sccs[0], m)
 
    logs = []
    for dis1 in merged_f.int_format():
        for dis2 in merged_f.int_format():
            if dis1 is not dis2 and list(set(dis1).intersection(set(dis2))): 
                for i in range(1, len(nempty_accs)):
                    logs.append(resolve_depend(aut, merged_f, nempty_sccs[0], nempty_accs[i], nempty_sccs[i]))                                                
    
    for i in range(1, len(nempty_accs)):
        l = (merge(aut, merged_f, nempty_sccs[0], nempty_accs[i], nempty_sccs[i]))   
        if all(log[0].states() != l[0].states() for log in logs):
            logs.append(l)
        else:
            for log in logs:
                if log[0] == l[0]:
                    log[1].update(l[1])
    return merged_f, logs 


### MAKE NEW AUT EQUIVALENT ###

def make_false(aut, scc, merged_acc):
    scc_clean_up_edges(aut, PACC(""), scc)
    add_m = []
    for dis in merged_acc.formula:
        if all(con.type == MarkType.Fin for con in dis):
            add_m.append(dis[0].num)
    
    if not add_m:
        return
    for s in scc.states():
        for e in aut.out(s):
            if e.dst in scc.states():
                for m in add_m:
                    e.acc.set(m)

def make_true(aut, scc, merged_acc):
    scc_clean_up_edges(aut, PACC(""), scc)

    min_inf = len(merged_acc[0])
    cheap_dis = merged_acc[0]
    for dis in merged_acc.formula:
        if all(con.type == MarkType.Fin for con in dis):
            return
        count = 0
        for con in dis:
            if con.type == MarkType.Inf:
                count += 1
        if count < min_inf:
            min_inf = count
            cheap_dis = dis

    add_m = []
    for con in cheap_dis:
        if con.type == MarkType.Inf:
            add_m.append(con.num)

    if not add_m:
        return
    for s in scc.states():
        for e in aut.out(s):
            if e.dst in scc.states():
                for m in add_m:
                    e.acc.set(m)


def make_equiv(aut, accs, sccs, logs, merged_acc):
    for i in range(len(accs)):
        if accs[i] is merged_acc:
            continue
        if accs[i].sat is True: 
            make_true(aut, sccs[i], merged_acc)
        elif (accs[i].sat is False) or (not accs[i].formula):
            make_false(aut, sccs[i], merged_acc)            
        else:
            contained = []
            for dis in accs[i].formula:
                for con in dis:
                    if con not in contained:
                        contained.append(con)
            for l in logs: 
                if l[0].states() == sccs[i].states():
                    for m in contained:
                        m.num = l[1][m.num]   
                    break     
            add_m = []
            for dis in merged_acc:
                if all(con not in contained for con in dis): #whole disjunct to be made false                    
                    if all(con.type == MarkType.Fin for con in dis):
                        if all(con.num not in add_m for con in dis):
                            add_m.append(dis[0].num)
                else: # all atomic conditions not present in disjunct to be made true
                    for con in dis:
                        if con not in contained and con.type == MarkType.Inf and con.num not in add_m:
                            add_m.append(con.num)
            for s in sccs[i].states():
                for e in aut.out(s):
                    if e.dst in sccs[i].states():
                        for m in add_m:
                            e.acc.set(m)
        scc_clean_up_edges(aut, merged_acc, sccs[i])


### MAIN ###

def eval_set(aut, mark, scc, m_all_e):
    if mark.num in m_all_e:
        return mark.type == MarkType.Inf
    if mark.num not in scc_current_marks(aut, scc):    
        return mark.type == MarkType.Fin
    return None        

def try_eval(aut, acc, scc):
    m_all_edges = set(scc_everywhere(aut, scc))

    eval_f = False
    for dis in acc.formula:
        eval_dis = True
        for con in dis:
            val = eval_set(aut, con, scc, m_all_edges)
            if eval_dis is False or val is False:
                eval_dis = False
            elif (eval_dis is None or eval_dis is True) and val is None:
                eval_dis = None
        if eval_dis is True:
            acc.sat = True
            acc.formula = []
            return
        elif (eval_f is None or eval_f is False) and eval_dis is None: 
            eval_f = None
    acc.set_sat(eval_f)
    if acc.sat is not None:
        acc.formula = []


def try_eval2(aut, acc, scc, weak):    
    if weak:
        if scc.is_rejecting():
            acc.set_sat(False)
        else:
            acc.set_sat(True)
        acc.formula = []
        return
    else:
        try_eval(aut, acc, scc)


def print_aut(aut, output, m):
    if output is not None:
        f = open(output, m)
        f.write(aut.to_str() + '\n')
        f.close()
    else:
        print(aut.to_str())


def process_aut(aut):
    spot.cleanup_acceptance_here(aut)
    if aut.get_acceptance().used_sets().count() < 1 or aut.prop_state_acc() == spot.trival.yes_value: 
        return

    accs = []
    sccs = []
    q = 0
    weak = spot.scc_info(aut).weak_sccs()
    for scc in spot.scc_info(aut):
        sccs.append(scc)
        acc = PACC(aut.get_acceptance().to_dnf())

        try_eval2(aut, acc, scc, weak[q])
        q += 1
        if acc.sat is None:
            simplify(aut, acc, scc)
        accs.append(acc)    
    
    new_acc, logs = merge_accs(aut, sccs, accs)
    if new_acc is None:
        if all(acc.sat is True for acc in accs):
            for scc in sccs:
                scc_clean_up_edges(aut, PACC(""), scc)
            aut.set_acceptance(0, spot.acc_code.t())
            return
        elif all(acc.sat is False for acc in accs):            
            for scc in sccs:
                scc_clean_up_edges(aut, PACC(""), scc)
            aut.set_acceptance(0, spot.acc_code.f())
            return
        else:
            new_acc = PACC("Inf(0)")
            for i in range(len(accs)):
                if accs[i].sat is True:
                    make_true(aut, sccs[i], new_acc)
                else:
                    make_false(aut, sccs[i], new_acc)
            aut.set_acceptance(1, spot.acc_code(str(new_acc)))
            return
    aut.set_acceptance(new_acc.max() + 1, spot.acc_code(str(new_acc))) 
    
    make_equiv(aut, accs, sccs, logs, new_acc)
    spot.cleanup_acceptance_here(aut)


def main(argv):    
    parser = argparse.ArgumentParser()
    parser.add_argument("-F", "--autfile", help="File containing automata in HOA format.")
    parser.add_argument("-O", "--outfile", help="File to print output to.")

    args = parser.parse_args()

    if not args.autfile: 
        print("No automata to process.", file=sys.stderr)

    aut = spot.automata(args.autfile)  
    for a in aut: 
        origin = spot.automaton(a.to_str())    
        process_aut(a) 

        if origin.get_acceptance().used_sets().count() < a.get_acceptance().used_sets().count():   
            a = origin            

        if args.outfile:
            print_aut(a, args.outfile, "a")
        else:
            print_aut(a, None, " ")               


if __name__ == "__main__":
   main(sys.argv[1:])
