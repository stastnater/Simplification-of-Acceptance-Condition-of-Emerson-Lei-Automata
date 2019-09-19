import spot
import enum
spot.setup()

### ACC CLASS ###
class MarkType(enum.Enum):
    Inf = 1
    Fin = 2

class ACCMark:
    def __init__(self, mtype, num):
        self.type = mtype
        self.num = num

    def __str__(self):
        return ''.join(["Inf" if self.type == MarkType.Inf else "Fin", '(', str(self.num), ')'])


class PACC:
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


### PARSE ACC ###
def parse_acc(acc):
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


### SIMPLIFY ###

def count_occur(acc, m): #count how many times mark m occurs in acc
    occur = 0
    for c in str(acc):
        if c == str(m):
            occur += 1
    return occur


def scc_current_marks(aut, scc): #return a list of marks currently present on edges in given scc
    marks = []
    for s in scc.states():
        for e in aut.out(s): 
            for m in e.acc.sets():
                if m not in marks:
                    marks.append(m)
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

       
def simplify_compl(aut, scc):
    for cm in scc_compl_marks(aut, scc):
        pass #TODO

### MERGE ACCs ###

### MAIN ###


### PARSE TESTS ###
def test_parser(filename):
    for aut in spot.automata(filename):
        print(aut.get_acceptance().to_dnf())
        new_formula = PACC(aut.get_acceptance().to_dnf())
        print("NEW: ", new_formula, "<-- ")



### RUN TESTS ###
FILENAME = '/home/tereza/Desktop/bp/allTela.aut'
test_parser(FILENAME)
