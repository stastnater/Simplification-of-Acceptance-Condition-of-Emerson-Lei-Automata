import spot
import enum
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


class PACC: #TODO: docu
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

    def int_format(self):
        f = []
        for dis in self.formula:
            d = []
            for con in dis:
                d.append(con.num)
            f.append(d)
        return f

    def clean_up(self, aut, scc):
        clean_f = []
        marks = scc_current_marks(aut, scc)
        for dis in self.formula:
            clean_dis = []
            for con in dis:
                if con.num in marks:
                    clean_dis.append(con)
            if not clean_dis:
                clean_f.append(clean_dis)
        self.formula = clean_f

    def get_mtype(self, m):
        for dis in self.formula:
            for con in dis:
                if m == con.num:
                    return con.type

    def find_m(self, m): #TODO: delete
        occurrences = []
        for i in range(len(self.formula) - 1):
            for j in range(len(self.formula[i]) - 1):
                if self.formula[i][j].num == m:
                    occurrences.append(i,j)
        return occurrences

    def find_m_dis(self, m):
        occurrences = []
        for i in range(len(self.formula) - 1):
            for j in range(len(self.formula[i]) - 1):
                if self.formula[i][j].num == m:
                    occurrences.append(i)
        return occurrences


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

def acc_count_occur(acc, m): 
    """
    Counts the occurrences of mark m in given acceptance condition.

    Parameters
    ----------
    acc : PACC
    m : int
        Number of an acceptance mark

    Returns
    -------
    int
        Amount of occurrences of m in acc.
    """

    occur = 0
    for dis in acc:
        for con in dis:
            if con.num == m:
                occur += 1
    return occur


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


def replace_marks(aut, scc, nm, rm): #replace rm (can contain multiple marks) with nm mark on all edges in scc
    for s in scc.states():
        for e in aut.out(s): 
            if e.dst in scc.states():
                for m in e.acc.sets():
                    if rm.has(m):
                        e.acc.clear(m)                    
                        e.acc.set(nm)


def scc_compl_marks(aut, scc): #return array of tuples of complementary marks in given scc
    c_marks = []
    for m1 in scc_current_marks(aut, scc):
        for m2 in scc_current_marks(aut, scc):
            if m1 is not m2:            
                are_compl = True
                for s in scc.states():
                    for e in aut.out(s):
                        if e.dst in scc.states and (m1 in e.acc.sets() and m2 in e.acc.sets()) or (m1 not in e.acc.sets() and m2 not in e.acc.sets()):
                            are_compl = False
                if are_compl and (m1, m2) not in c_marks and (m2, m1) not in c_marks:
                    c_marks.append((m1, m2))
    return c_marks


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
    

def simpl_inf_con(autacc, scc, subsets):
    int_acc = autacc[1].int_format()
    for sub in subsets:
        if autacc[1].get_mtype(sub[0]) == 1 and autacc[1].get_mtype(sub[1]) == 1:
            if autacc[1].find_m_dis(sub[0]) == autacc[1].find_m_dis(sub[1]):
            #TODO: remove all sub[0] marks from aut, clean acc, ret autacc
                pass
            else:
                pass #remove sub[0] at least from dis where it is with sub[1]
        
                



def simpl_subsets(aut, scc, acc):
    subsets = scc_subsets(aut, scc)
    autacc = (aut, acc)
       
def simplify_compl(aut, scc):
    for cm in scc_compl_marks(aut, scc):
        pass #TODO:



### SIMPLIFY ###

def simplify(aut, scc, acc):
    
    pass

### MERGE ACCs ###

### MAIN ###

def main():
    aut = spot.automaton('necofile.aut')
    scc_accs = [aut.get_acceptance().to_dnf()]*(spot.scc_info(aut).scc_count())
    acc_index = 0

    for scc in spot.scc_info(aut):
        simplify(aut, scc, acc_index)
        acc_index += 1


### PARSE TESTS ###
def test_parser(filename):
    for aut in spot.automata(filename):
        print(aut.get_acceptance().to_dnf())
        new_formula = PACC(aut.get_acceptance().to_dnf())
        print("NEW: ", new_formula, "<-- ")



### RUN TESTS ###
FILENAME = '/home/tereza/Desktop/bp/allTela.aut'
test_parser(FILENAME)
