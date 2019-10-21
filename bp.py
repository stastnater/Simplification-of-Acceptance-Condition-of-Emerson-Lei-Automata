import spot
import enum
import sys
import copy
from scipy.optimize import linear_sum_assignment
spot.setup()

### ACC CLASS ###
class MarkType(enum.Enum): #TODO: docu
    Inf = 1
    Fin = 2


class ACCMark: #TODO: docu
    def __init__(self, mtype, num):
        self.type = mtype
        self.num = num

    def __str__(self):
        return ''.join(["Inf" if self.type == MarkType.Inf else "Fin", '(', str(self.num), ')'])

    def __eq__(self, other): 
        if (self.type == other.type) and (self.num == other.num): 
            return True
        else: 
            return False


class PACC: #TODO: docu #NO OH GOD NO NO NOOOOOOOOOO
    def __init__(self, acc):
        self.formula = parse_acc(acc)

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
        unit_dis = []
        for dis in self.formula:
            if len(dis) == 1:
                unit_dis.append(dis[0])

        res_f = []
        for dis in self.formula:
            if len(dis) == 1:
                res_f.append(dis)
            else:
                if all(con not in unit_dis for con in dis):
                    res_f.append(dis)
                    print("found redundant dis")
                   
        self.formula = []
        for dis in res_f:
            if dis not in self.formula:
                self.formula.append(dis)

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
        #print("found ", m, " at ", occurrences)              
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
    """
    Parses an acc in DNF and returns the formula represented by list of lists of ACCMarks. 
    The inner lists represent disjuncts of the formula. 
    The inner lists contain ACCMark (see ACCMark class documentation) objects representing atomic conditions (such as Inf(1)).

    Example: (Fin(1) & Inf(2)) | (Inf(3)) | (Fin(1) & Fin(4)) --> [[ACCMark(2,1), ACCMark(1,2)] [ACCMark(1,3)] [ACCMark(2,1), ACCMark(2,4)]]

    Parameters
    ----------
    aut : spot::acc_cond::acc_code

    Returns
    -------
    [[ACCMark]]
        List of lists of ACCMarks.
    """

    formula = []
    for dis in str(acc).split('|'):
        new_dis = []
        for con in dis.split('&'):
            mtype = None
            if con.find("Fin") != -1:
                mtype = MarkType.Fin
            else:
                mtype = MarkType.Inf
            for c in con:
                if c.isdigit():                    
                    new_dis.append(ACCMark(mtype, int(c)))
        formula.append(new_dis)            
    return formula


### SIMPLIFY AUXILIARY ###

def scc_clean_up_edges(aut, acc, scc):
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
    """
    Return a list of marks currently present on edges in given scc.

    Parameters
    ----------
    aut : spot::twa        
    scc : spot::scc_info_node 

    Returns
    -------
    [int]
        List of marks on edges of given scc.
    """

    marks = []
    for s in scc.states():
        for e in aut.out(s): 
            for m in e.acc.sets():
                if m not in marks:
                    marks.append(int(m)) 
    return marks


def scc_compl_sets(aut, scc): #return array of tuples of complementary marks in given scc #TODO: docu
    compl = []
    marks = scc_current_marks(aut, scc)
    for m1 in marks:
        for m2 in marks:
            if m1 is not m2:            
                are_compl = True
                for s in scc.states():
                    for e in aut.out(s):
                        if e.dst in scc.states() and (m1 in e.acc.sets() and m2 in e.acc.sets()) or (m1 not in e.acc.sets() and m2 not in e.acc.sets()):
                            are_compl = False
                if are_compl and (m1, m2) not in compl and (m2, m1) not in compl:
                    compl.append((m1, m2))
    return compl



def scc_subsets(aut, scc):
    """
    Return a list of tuples (m1, m2) containing acceptance marks m1 and m2 where m2 is a subset of m1 in given scc.

    Parameters
    ----------
    aut : spot::twa        
    scc : spot::scc_info_node 

    Returns
    -------
    [(int, int)]
        List of tuples of subsets.
    """

    marks = scc_current_marks(aut, scc)
    subsets = []
    for m1 in marks:
        for m2 in marks:
            if m1 is not m2:
                is_sub = True
                for s in scc.states():
                    for e in aut.out(s):
                        if e.dst in scc.states() and e.acc.has(m2) and not e.acc.has(m1):
                            is_sub = False
                if is_sub:
                    subsets.append((m1, m2))
    return subsets 
 


def remove_mark(aut, scc, m): 
    """
    Removes mark m from all edges in the given scc.

    Parameters
    ----------
    aut : spot::twa        
    scc : spot::scc_info_node 
    m   : int
    """

    for s in scc.states():
        for e in aut.out(s): 
            if e.dst in scc.states():
                if e.acc.has(m):
                    e.acc.clear(m)
                     
                        
def simpl_inf_con(aut, acc, scc, subsets): #TODO: docu
    for sub in subsets:        
        if acc.get_mtype(sub[0]) == MarkType.Inf and acc.get_mtype(sub[1]) == MarkType.Inf:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            if (all(i in sub_i for i in super_i)):
                #print("subset in: ", sub_i, "   superset in: ", super_i)
                print("inf con removing: ", sub[0])
                remove_mark(aut, scc, sub[0])
                acc.clean_up(aut, scc)                            
            else:
                for i in super_i:
                    if i in sub_i:
                        acc.rem_from_dis(i, sub[0])
                        acc.clean_up(aut, scc)


def simpl_fin_con(aut, acc, scc, subsets):
    #TODO: count dis with fin, if discount < unique fin -> merge all fins in each dis
    # else: find all pairs which always appear together -> merge them -> recursively?
    for sub in subsets:
        if acc.get_mtype(sub[0]) == MarkType.Fin and acc.get_mtype(sub[1]) == MarkType.Fin:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            if (all(i in super_i for i in sub_i)):
                print("fin con removing: ", sub[1])
                remove_mark(aut, scc, sub[1])
                acc.clean_up(aut, scc)
            else:
                for i in sub_i:
                    if i in super_i:
                        acc.rem_from_dis(i, sub[1])
                        acc.clean_up(aut, scc)
        
                
def simpl_inf_dis(aut, acc, scc):
    unit_inf = []
    for dis in acc.formula:
        if len(dis) == 1 and dis[0].type == MarkType.Inf:
            unit_inf.append(dis)
    if len(unit_inf) < 2:
        return
    for i in range(1, len(unit_inf) - 1):
        replace_marks(aut, scc, unit_inf[i][0].num, unit_inf[0][0].num)
    acc.clean_up(aut, scc)            

    """ if restored add subsets to param
    for sub in subsets:
        if acc.get_mtype(sub[0]) == MarkType.Inf and acc.get_mtype(sub[1]) == MarkType.Inf:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            if (any(len(acc[i]) == 1 for i in sub_i)) or (any(len(acc[i]) == 1 for i in super_i)):
                print("inf dis removing: ", sub[1])
                remove_mark(aut, scc, sub[1])
                acc.clean_up(aut, scc)
    """

def simpl_fin_dis(aut, acc, scc, subsets):
    for sub in subsets:
        if acc.get_mtype(sub[0]) == MarkType.Fin and acc.get_mtype(sub[1]) == MarkType.Fin:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            if (any(len(acc[i]) == 1 for i in sub_i)) or (any(len(acc[i]) == 1 for i in super_i)):
                print("fin dis removing: ", sub[0])
                remove_mark(aut, scc, sub[0])
                acc.clean_up(aut, scc)


def simpl_co_dis(aut, acc, scc, compl_sets): 
    for co in compl_sets:
        if acc.get_mtype(co[0]) != acc.get_mtype(co[1]):
            inf, fin = co[0], co[1] 
            if acc.get_mtype(co[0]) == MarkType.Fin:
                inf, fin = co[1], co[0]
            inf_i = acc.find_m_dis(inf)
            fin_i = acc.find_m_dis(fin)
            if all(i not in fin_i for i in inf_i):
                remove_mark(aut, scc, fin)
                acc.clean_up(aut, scc)


def simpl_co_con(aut, acc, scc, compl_sets): #TODO: is redundant?
    for co in compl_sets:
        if acc.get_mtype(co[0]) != acc.get_mtype(co[1]):
            inf, fin = co[0], co[1] 
            if (acc.get_mtype(co[0]) == MarkType.Fin):
                 inf, fin = co[1], co[0]
            inf_i = acc.find_m_dis(inf)
            fin_i = acc.find_m_dis(fin)
            if all(i in fin_i for i in inf_i):
                remove_mark(aut, scc, inf)
                acc.clean_up(aut, scc)
            else:
                for i in inf_i:
                    if i in fin_i:
                        acc.rem_from_dis(i, inf)
                        acc.clean_up(aut, scc)


def simpl_false(aut, acc, subsets, compl_sets):
    for sub in subsets:
        if acc.get_mtype(sub[0]) == MarkType.Fin and acc.get_mtype(sub[1]) == MarkType.Inf:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            for i in sub_i:
                if i in super_i:
                    acc.rem_dis(i)
                    print("false dis removed at ", i)
    for co in compl_sets:
        if acc.get_mtype(co[0]) == MarkType.Fin and acc.get_mtype(co[1]) == MarkType.Fin:
            co1_i = acc.find_m_dis(co[0])
            co2_i = acc.find_m_dis(co[1])
            for i in co1_i:
                if i in co2_i:
                    acc.rem_dis(i)
                    print("false dis removed at ", i)


### SIMPLIFY ###

def simplify(aut, acc, scc):
    acc_l = acc.acc_len()
    acc.clean_up(aut, scc)
    subsets = scc_subsets(aut, scc)
    compl_sets = scc_compl_sets(aut, scc)
    print("SCC: ", scc.states())
    print("starting acc: ", acc)

    simpl_inf_con(aut, acc, scc, subsets)
    simpl_fin_con(aut, acc, scc, subsets)
    simpl_inf_dis(aut, acc, scc)
    simpl_fin_dis(aut, acc, scc, subsets)
    simpl_false(aut, acc, subsets, compl_sets)
    simpl_co_con(aut, acc, scc, compl_sets)
    simpl_co_dis(aut, acc, scc, compl_sets)
    scc_clean_up_edges(aut, acc, scc)
    print(acc)

    if acc_l > acc.acc_len():
        simplify(aut, acc, scc)
    

### MERGE AUXILIARY ###

def shift_fst_acc(aut, acc, scc):
    next_m = aut.acc().num_sets()
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
    return log


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
            if e.dst in scc.states() and e.acc.has(origin_m):
                e.acc.clear(origin_m)
                if not e.acc.has(new_m):
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


def place_con(dis, con, used):
    if con in dis and not used[dis.index(con)]:
        return dis.index(con)
    for i in range(len(dis)):
        if dis[i].type == con.type and not used[i]:            
            return i
    return None
        

def add_dupl_marks(aut, scc, origin_m, new_m):
    for s in scc.states():
        for e in aut.out(s):
            if e.dst in scc.states() and e.acc.has(origin_m):
                if not e.acc.has(new_m):
                    e.acc.set(new_m)


def get_dependencies(acc):
    depend = []
    int_f = acc.int_format()
    for i in range(len(int_f)):        
        for j in range(len(int_f)):
            intersect = list(set(int_f[i]).intersection(set(int_f[j])))
            if i != j and intersect:
                for mark in intersect:
                    if any(d[0].num == mark for d in depend):
                        for d in depend:
                            if d[0].num == mark:
                                if i not in d[1:]:
                                    d.append(i)
                                if j not in d[1:]:
                                    d.append(j)                                    
                    else:
                        depend.append([ACCMark(acc.get_mtype(mark), mark), i, j])
    return depend

def resolve_dependencies(aut, acc1, scc1, acc2, scc2):
    acc1_d = get_dependencies(acc1)
    m = make_matrix(acc1, acc2)
    row_ind, col_ind = linear_sum_assignment(m) 
    s_acc2 = []
    for ci in col_ind:
        if ci < len(acc2):
            s_acc2.append(acc2[ci])
        else: 
            s_acc2.append([])
    acc2_d = get_dependencies(PACC(str(s_acc2)))
    for d1 in acc1_d: #TODO: sort acc1_d by length for better results
        for d2 in acc2_d:
            # dependencies needn't be removed
            if d1[0].type == d2[0].type and all(d1[i] in d2 for i in range(1, len(d1))):
                add_dupl_marks(aut, scc2, d2[0].num, d1[0].num) #TODO: possibly replace, not duplicate
                """
                for s in scc2.states():
                    for e in aut.out(s):
                        if e.acc.has(d2[0].num):
                            e.acc.set(d1[0].num) 
                """
                for dis in acc2:
                    for con in dis:
                        if con == d2[0]:
                            con.num = d1[0].num           
                d2 = []
                d1 = []
                break 
    for d1 in acc1_d:
        if d1:
            for i in range(2, len(d1)):
                for con in acc1[d1[i]]:
                    if con == d1[0]:
                        new_m = acc1.max() + 1
                        add_dupl_marks(aut, scc1, con.num, new_m)
                        con.num = new_m
                                                         

### MERGE ###

def merge(aut, acc1, scc1, acc2, scc2):
    m = make_matrix(acc1, acc2)
    print("accs: ", acc1, " and ", acc2)
    print(m)
    row_ind, col_ind = linear_sum_assignment(m) 
    print("row ind: ", row_ind, " col ind: ", col_ind)

    log = [scc2]
    for i in range(len(row_ind)):
        dis1 = acc1[col_ind[i]]
        used = [False] * len(dis1)
        dis2 = []
        if row_ind[i] < len(acc2):
            dis2 = acc2[row_ind[i]]
        print("will merge: ", row_ind[i], col_ind[i])
        for con in dis2:
            index = place_con(dis1, con, used)
            if index is None:
                log.append((con.num, acc1.max() + 1))
                dis1.append(ACCMark(con.type, acc1.max() + 1))
                used.append(True)
            else:
                log.append((con.num, dis1[index].num))
                replace_marks(aut, scc2, con.num, dis1[index].num)
                used[index] = True
    return log
                        
def merge_accs(aut, sccs, accs):
    nempty_accs = []
    nempty_sccs = []

    for i in range(len(accs)):
        if accs[i].formula:
            nempty_accs.append(accs[i])
            nempty_sccs.append(sccs[i])

    acc_quicksort(nempty_accs, nempty_sccs, 0, len(nempty_accs) - 1)
    
    merged_f = PACC(str(nempty_accs[0]))
    merged_log = shift_fst_acc(aut, merged_f, nempty_sccs[0])
    
    for dis1 in merged_f.int_format():
        for dis2 in merged_f.int_format():
            if dis1 != dis2 and list(set(dis1).intersection(set(dis2))):
                for i in range(1, len(nempty_accs)):
                    resolve_dependencies(aut, merged_f, nempty_sccs[0], nempty_accs[i], nempty_sccs[i])
                break
    
    logs = []
    for i in range(1, len(nempty_accs) - 1):
        logs.append(merge(aut, merged_f, nempty_sccs[0], nempty_accs[i], nempty_sccs[i]))
    return merged_f, merged_log, logs


### MAKE NEW AUT EQUIVALENT ###

def make_false(aut, scc, merged_acc):
    for s in scc.states(): #TODO: clean up instead?
        for e in aut.out(s):
            if e.dst in scc.states():
                for m in e.acc.sets():
                    e.acc.clear(m)
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

def make_equiv(aut, accs, sccs, log, logs, merged_acc):
    print(log)
    for i in range(1, len(accs)):
        if not accs[i]: 
            make_false(aut, sccs[i], merged_acc)
        else:
            contained = []
            for dis in accs[i].formula:
                for con in dis:
                    if con not in contained:
                        contained.append(con)
            for l in logs: 
                if l[0] is sccs[i]:
                    for m in contained:
                        m.num = l[m.num]   
                    break     
            add_m = []
            for dis in merged_acc:
                if all(con not in contained for con in dis): #whole disjunct to be made false
                    if all(con.type == MarkType.Fin for con in dis):
                        if all(con.num not in add_m for con in dis):
                            add_m.append(dis[0].num)
                else: # all atomic conditions not present in disjunct to be made true
                    for con in dis:
                        if con not in contained and con.type == MarkType.Inf:
                            add_m.append(con.num)
            for s in sccs[i].states():
                for e in aut.out(s):
                    if e.dst in sccs[i].states():
                        for m in add_m:
                            e.acc.set(m)
        scc_clean_up_edges(aut, merged_acc, sccs[i])

### MAIN ###

def main(argv):
    FILENAME = str(sys.argv[1])
    aut = spot.automaton("/home/tereza/Desktop/bp/" + FILENAME)
    accs = []
    sccs = []

    #print(aut.get_acceptance().to_dnf())
    for scc in spot.scc_info(aut):
        sccs.append(scc)
        acc = PACC(aut.get_acceptance().to_dnf())

        simplify(aut, acc, scc)
        accs.append(acc)

    new_acc, log, logs = merge_accs(aut, sccs, accs)
    aut.set_acceptance(new_acc.max() + 1, spot.acc_code(str(new_acc)))

    make_equiv(aut, accs, sccs, log, logs, new_acc)
    
    aut.save('_' + FILENAME)
    print(new_acc)
    


"""
### PARSE TESTS ###
def test_parser(filename):
    for aut in spot.automata(filename):
        print(aut.get_acceptance().to_dnf())
        new_formula = PACC(aut.get_acceptance().to_dnf())
        print("NEW: ", new_formula, "<-- ")



### RUN TESTS ###
FILENAME = '/home/tereza/Desktop/bp/allTela.aut'
test_parser(FILENAME)
"""

if __name__ == "__main__":
   main(sys.argv[1:])
