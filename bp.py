import spot
import enum
import sys
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
                
    def int_format(self):
        f = []
        for dis in self.formula:
            d = []
            for con in dis:
                d.append(con.num)
            f.append(d)
        return f

    def resolve_redundancy(self):
        """
        int_f = self.int_format()
        alone_m = []
        for dis in int_f:
            if len(dis) == 1:
                alone_m.append(dis[0])
        new_f = []
        for i in range(len(int_f) - 1):
            if len(int_f[i]) == 1 or all(con not in alone_m for con in int_f[i]):
                new_f.append(self.formula[i])        
        
        self.formula = []
        for dis in new_f:
            if dis not in self.formula:
                self.formula.append(dis)
                """
        res_f = self.formula
        for dis in self.formula:
            if len(dis) == 1:
                for dis2 in res_f:
                    if dis[0] in dis2:
                        dis2 = dis
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
        print("found ", m, " at ", occurrences)              
        return occurrences

    def rem_from_dis(self, index, m):
        new_dis = []
        for con in self.formula[index]:
            if con.num != m:
                new_dis.append(con)
        self.formula[index] = new_dis

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
                print("subset in: ", sub_i, "   superset in: ", super_i)
                print("inf con removing: ", sub[0])
                remove_mark(aut, scc, sub[0])
                acc.clean_up(aut, scc)                            
            else:
                for i in super_i:
                    if i in sub_i:
                        acc.rem_from_dis(i, sub[0])


def simpl_fin_con(aut, acc, scc, subsets): #TODO: docu
    for sub in subsets:
        if acc.get_mtype(sub[0]) == MarkType.Fin and acc.get_mtype(sub[1]) == MarkType.Fin:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            if (all(i in sub_i for i in super_i)):
                print("fin con removing: ", sub[1])
                remove_mark(aut, scc, sub[1])
                acc.clean_up(aut, scc)
            else:
                for i in sub_i:
                    if i in super_i:
                        acc.rem_from_dis(i, sub[1])
        
                
def simpl_inf_dis(aut, acc, scc, subsets):
    for sub in subsets:
        if acc.get_mtype(sub[0]) == MarkType.Inf and acc.get_mtype(sub[1]) == MarkType.Inf:
            sub_i = acc.find_m_dis(sub[1])
            super_i = acc.find_m_dis(sub[0])
            if (all(len(acc[i]) == 1 for i in sub_i)) and (all(len(acc[i]) == 1 for i in super_i)):
                print("inf dis removing: ", sub[1])
                remove_mark(aut, scc, sub[1])
                acc.clean_up(aut, scc)


def simplify_compl(aut, scc):
    for cm in scc_compl_marks(aut, scc):
        pass #TODO: ???


### SIMPLIFY ###

def simplify(aut, scc, acc):
    acc_l = acc.acc_len()
    subsets = scc_subsets(aut, scc)
    acc.clean_up(aut, scc)
    print("SCC: ", scc.states())
    print("starting acc: ", acc, "     int format: ", acc.int_format())

    simpl_inf_con(aut, acc, scc, subsets)
    print(acc)
    simpl_fin_con(aut, acc, scc, subsets)
    print(acc)
    simpl_inf_dis(aut, acc, scc, subsets)
    print(acc)



    if acc_l > acc.acc_len():
        #simplify(aut, scc, acc) TODO: uncomment when done, or iterative?
        pass

    

### MERGE ACCs ###

### MAIN ###

def main(argv):
    FILENAME = str(sys.argv[1])
    aut = spot.automaton("/home/tereza/Desktop/bp/" + FILENAME)
    scc_accs = []

    #print(aut.get_acceptance().to_dnf())

    for scc in spot.scc_info(aut):
        acc = PACC(aut.get_acceptance().to_dnf())

        print(acc, "\nlen: ", len(acc))
        print(acc.int_format)
        print("index 0: ", acc[0])

        simplify(aut, scc, acc)
        scc_accs.append(acc)
    
    aut.save('_' + FILENAME)

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
