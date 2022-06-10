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

myYield = ROOT.RooRealVar("yield","yield",50,0,100)
myObs = ROOT.RooRealVar("x","x",5,0,10)
myPar1 = ROOT.RooRealVar("y","y",5,3,6)
myPar2 = ROOT.RooRealVar("z","z",1,0,5)
myGaus = ROOT.RooGaussian("myPdf","myPdf",myObs,myPar1,myPar2)
myModel = ROOT.RooExtendPdf("myModel","myModel",myGaus,myYield)

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
    if objName=="myData":
        w = ROOT.RooRealVar("weightVar","weight",1)
        c = ROOT.RooCategory("channelCat","channelCat")
        c.defineType("SR");c.defineType("CR")
        s = ROOT.RooArgSet(mySetVars[0],w,c);
        d = ROOT.RooDataSet("myData","My Dataset",s,"weightVar")
        mySetVars[0].setVal(3)
        d.add(s,30)
        c.setIndex(1)
        d.add(s,20)
        mySetVars[0].setVal(2)
        d.add(s,15)
        g1 = ROOT.RooRealVar("glob1","glob1",0); g2 = ROOT.RooRealVar("another_glob","another global obs",2)
        d.setGlobalObservables(ROOT.RooArgSet(g1,g2))
        return d
    if objName=="myModel":
        return myModel
    if objName=="myWorkspace":
        ROOT.RooMsgService.instance().getStream(ROOT.RooFit.INFO).removeTopic(ROOT.RooFit.NumIntegration)
        def own(obj):
            ROOT.SetOwnership(obj,True)
            return obj
        f = ROOT.TFile("tutorialModelTemplate.root")
        w = f.Get("combined")
        model = w.pdf("simPdf")
        obs = own(own(own(model.getVariables()).selectByAttrib("obs",True)).selectByAttrib("global",False))
        globs = own(own(model.getVariables()).selectByAttrib("global",True))
        w.var("sig_mass").setVal(100)
        w.var("mu").setVal(0)
        d = model.generate(obs, ROOT.RooFit.Extended())
        d.SetName("obsData")
        globs.first().setVal(5.5)
        d.setGlobalObservables(globs)
        w.Import(d,ROOT.RooFit.Silence(True))
        #w.saveSnapshot("obsData",globs)
        globs.first().setVal(5)
        return w
    
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
    
    
def test_3a(ans):
    try:
        if ans != "channelCat": printRed("test_3a: INCORRECT - wrong variable name")
        else: printGreen("test_3a: CORRECT")
    except:
        printRed("test_3a: INCORRECT - didn't give a string")
def test_3b(ans):
    try:
        if ans != 15: printRed("test_3b: INCORRECT - wrong weight")
        else: printGreen("test_3b: CORRECT")
    except:
        printRed("test_3b: INCORRECT - didn't give a number")       
def test_3c(ans):
    try:
        if ans != 2: printRed("test_3c: INCORRECT - wrong number of global observables")
        else: printGreen("test_3c: CORRECT")
    except:
        printRed("test_3c: INCORRECT - didn't give a number")  
        
def test_4a(ans):
    try:
        def own(obj):
            ROOT.SetOwnership(obj,True)
            return obj
        graph = ROOT.TGraph()
        var = own(myModel.getVariables()).find("x")
        import numpy as np
        for x in np.arange(var.getMin(),var.getMax(),0.2):
            var.setVal( x )
            graph.AddPoint(var.getVal(),myModel.getVal(ROOT.RooArgSet(var)))
        if ans.GetN()!=graph.GetN(): 
            printRed("test_4a: INCORRECT - wrong number of points")
            return
        
        for i in range(graph.GetN()):
            if (abs(graph.GetPointX(i)-ans.GetPointX(i))>0.0001):
                printRed("test_4a: INCORRECT - wrong x value: {}".format(ans.GetPointX(i)))
                return
            elif (abs(graph.GetPointY(i)-ans.GetPointY(i))>0.0001):
                printRed("test_4a: INCORRECT - wrong y value: {}".format(ans.GetPointY(i)))
                return
        
        else: printGreen("test_4a: CORRECT")
    except:
        printRed("test_4a: INCORRECT - didn't give a graph")
        
def test_5a(ans):
    try:
        if ans != "simPdf": printRed("test_5a: INCORRECT - wrong pdf name")
        else: printGreen("test_5a: CORRECT")
    except:
        printRed("test_5a: INCORRECT - didn't give a string") 
def test_5b(ans):
    try:
        if ans != "obsData": printRed("test_5b: INCORRECT - wrong dataset name")
        else: printGreen("test_5b: CORRECT")
    except:
        printRed("test_5b: INCORRECT - didn't give a string") 
def test_5c(name,val):
    try:
        if name != "globs_alpha_par": printRed("test_5c: INCORRECT - wrong global observable name")
        elif val != 5.5: printRed("test_5c: INCORRECT - wrong value")
        else: printGreen("test_5c: CORRECT")
    except:
        printRed("test_5c: INCORRECT - didn't give a name and number")         