

class Clause:
    def __init__(self, pred1, pred2,neg1,neg2):
        self.pred1=pred1
        self.pred2=pred2
        self.neg1=neg1
        self.neg2=neg2

    def __eq__(self, obj):
        if type(self)!=type(obj):
            return False

        return self.pred1==obj.pred1 and self.pred2==obj.pred2 and self.neg2==obj.neg2 and self.neg1==obj.neg1

class Pred:
    def __init__(self, symbol, arg1,arg2):
        self.symbol=symbol
        self.arg1=arg1
        self.arg2=arg2

    def copy(self,p):
        self.symbol=p.symbol
        if type(p.arg1) == type("s"):
            self.arg1 = p.arg1
        elif type(p.arg1) == type(self):
            self.arg1.symbol = p.arg1.symbol
            self.arg1.arg1 = p.arg1.arg1
            self.arg1.arg2 = p.arg1.arg2

        self.arg2 = p.arg2

    def __eq__(self, obj):
        if type(self)!=type(obj):
            return False
        if type(self.arg1)!=type(obj.arg1) or type(self.arg2)!=type(obj.arg2):
            return False

        return self.symbol==obj.symbol and self.arg1==obj.arg1 and self.arg2==obj.arg2
def s_to_pred(s):
    symbol=s[0]
    if s[3]=="(":
        arg1=Pred(s[2],s[4],None)
    else:
        arg1=s[2]
    if "," in s:
        i=s.index(",")
        if s[i+2] == "(":
            arg2 = Pred(s[i+1], s[i+3], None)
        else:
            arg2 = s[i+1]
    else:
        arg2=None

    return Pred(symbol,arg1,arg2)


def pred_to_s(pred):
    s=pred.symbol+"("
    if type(pred.arg1)==type("s"):
        s+=pred.arg1
    else:
        s += pred.arg1.symbol+"("+pred.arg1.arg1+")"

    if pred.arg2 == None:
        s += ")"
    elif type(pred.arg2)==type("s"):
        s += "," + pred.arg2 + ")"
    elif type(pred.arg2)==type(pred):
        s += "," + pred.arg2.symbol + "(" + pred.arg2.arg1 + ")"+")"

    return s
def s_to_clause(s):
    neg1 = False
    neg2 = None
    if "+" in s:
        i=s.index("+")
        if s[0]=="~":
            neg1 = True
            pred1 = s_to_pred(s[1:i])
        else:
            pred1 = s_to_pred(s[0:i])
        if s[i+1]=="~":
            neg2 = True
            pred2 = s_to_pred(s[i+2:])
        else:
            neg2 = False
            pred2 = s_to_pred(s[i+1:])

        return Clause(pred1, pred2, neg1, neg2)

    else:
        if s[0]=="~":
            neg1 = True
            pred1 = s_to_pred(s[1:])
        else:
            pred1 = s_to_pred(s[0:])

        return Clause(pred1, None, neg1, neg2)

def clause_to_s(clause):
    s=""

    if(clause.neg1):
        s+="~"

    s += pred_to_s(clause.pred1)
    if clause.neg2!=None:
        s+="+"

        if (clause.neg2):
            s += "~"
        s += pred_to_s(clause.pred2)
    return s
def is_neg(clause,obt):
    if(obt=="empty"):
        return False
    if obt.pred1.symbol == clause.pred1.symbol:
        if obt.neg1 == True and clause.neg1 == False or obt.neg1 == False and clause.neg1 == True:
            return True


    elif clause.pred2!=None and obt.pred1.symbol == clause.pred2.symbol:
        if obt.neg1 == True and clause.neg2 == False or obt.neg1 == False and clause.neg2 == True:
            return True
    return False

def conclude(clause1,clause2):
    pred=Pred(None,None,None)
    neg=None
    if(clause2.neg2==None and clause1.neg2==None):
        return "empty"
    elif clause1.neg2 == None:

        if clause1.pred1.symbol==clause2.pred1.symbol:
            neg=clause2.neg2
            pred.copy(clause2.pred2)

        else:
            neg = clause2.neg1
            pred.copy(clause2.pred1)

    elif clause2.neg2 == None:
        if clause2.pred1.symbol == clause1.pred1.symbol:
            neg = clause1.neg2
            pred.copy(clause1.pred2)

        else:
            neg = clause1.neg1
            pred.copy(clause1.pred1)

    if clause1.pred1.arg1.isupper():
        pred.arg1=clause1.pred1.arg1
    if clause2.pred1.arg1.isupper():
        pred.arg1 = clause2.pred1.arg1
    return Clause(pred,None,neg,None)

def theorem_prover(base_clauses,obt_clauses):
    temp=[]

    for e in base_clauses:
        temp.append(s_to_clause(e))

    base_clauses=temp
    temp=[]

    for e in obt_clauses:
        temp.append(s_to_clause(e))

    obt_clauses=temp
    flag=True
    retList=[]
    while flag:
        flag=False
        temp=[]
        for e in obt_clauses:
            for d in base_clauses:
                if is_neg(d,e):
                    temp.append((e,d))

        for e in temp:
            c=conclude(e[0],e[1])
            if c not in obt_clauses:
                obt_clauses.append(conclude(e[0],e[1]))

            if c=="empty":
                flag=False
                s=clause_to_s(e[0]) + "$" + clause_to_s(e[1]) + "$" + c
                if s not in retList:
                    retList.append(s)
            else:
                s=clause_to_s(e[0]) + "$" + clause_to_s(e[1]) + "$" + clause_to_s(c)
                if s not in retList:
                    retList.append(s)
        if temp==[]:
            break

    if("empty" not in retList[-1] ):
        return("no", [])
    else:
        return("yes",retList)




