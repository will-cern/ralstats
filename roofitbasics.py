import sys
if "/usr/local/root/lib" not in sys.path: sys.path.append("/usr/local/root/lib")
import ROOT

ROOT.RooFit # call here just to get banner out of the way

def printRed(text):
    print('\033[91m' + text + '\033[0m')
def printGreen(text):
    print('\033[92m' + text + '\033[0m')

myVar = None
mySet = None
myList = None

# define here to keep alive
mySetVars = [ROOT.RooRealVar("foo","",3),ROOT.RooRealVar("bar","",6,0,10),ROOT.RooCategory("baz","")]
myListVars = [ROOT.RooRealVar("foo","",4),ROOT.RooCategory("baz","")]
    
def getObject(objName):
    global myVar,mySet,myList
    if objName=="myVar":
        myVar = ROOT.RooRealVar("myVar","my variable",5)
        return myVar
    if objName=="mySet":
        s = ROOT.RooArgSet(); s.setName("mySet")
        for a in mySetVars: s.add(a)
        mySet = s
        return s
    if objName=="myList":
        s = ROOT.RooArgList(); s.setName("myList")
        for a in myListVars: s.add(a)
        myList = s
        return s
    
    raise RuntimeError("Unknown object {}".format(objName))
    
def test_1a(ans):
    if ans==True: printGreen("test_1a: CORRECT")
    else: printRed("test_1a: INCORRECT")
        
def test_1b(r):
    try:
        if r.getVal()!=7: printRed("test_1b: INCORRECT - wrong value")
        elif r.getError()!=4: printRed("test_1b: INCORRECT - wrong error")
        elif r.getMin()!=0 or r.getMax()!=10: printRed("test_1b: INCORRECT - wrong range")
        else: printGreen("test_1b: CORRECT")
    except:
        printRed("test_1b: INCORRECT - didn't give a RooRealVar")
    
    
def test_2a(ans):
    try:
        if ans=="baz": printGreen("test_2a: CORRECT")
        else: printRed("test_2a: INCORRECT")
    except:
        printRed("test_2a: INCORRECT - didn't give a string")
def test_2b(ans):
    try:
        if ans==6: printGreen("test_2b: CORRECT")
        else: printRed("test_2b: INCORRECT")
    except:
        printRed("test_2b: INCORRECT - didn't give a number")
def test_2c(ans):
    try:
        if ans.size()!=2 or ans.find("bar")==None or ans.find("baz")==None: printRed("test_2c: INCORRECT")
        else: printGreen("test_2c: CORRECT")
    except:
        printRed("test_2c: INCORRECT - didn't give a collection")
def test_2d(ans):
    try:
        if ans.size()!=1 or ans.find("bar")==None: printRed("test_2d: INCORRECT")
        else: printGreen("test_2d: CORRECT")    
    except:
        printRed("test_2d: INCORRECT - didn't give a collection")