import sys
if "/usr/local/root/lib" not in sys.path: sys.path.append("/usr/local/root/lib")
import ROOT
import tqdm
import math

if not hasattr(ROOT,"Asymptotica"):
    ROOT.gInterpreter.ProcessLine("#include \"Asymptotics.h\"")
    ROOT.Asymptotica
    
ROOT.gInterpreter.ProcessLine(".X style.C")
ROOT.RooMsgService.instance().getStream(ROOT.RooFit.INFO).removeTopic(ROOT.RooFit.NumIntegration)

def own(obj):
    ROOT.SetOwnership(obj,True)
    return obj

def plot(model,dataset,globs,fitResult=None):
    oldPad = ROOT.gPad
    if not ROOT.gROOT.GetListOfCanvases().FindObject("plot"):
        c = ROOT.TCanvas("plot","",700,215)
        ROOT.SetOwnership(c,False)
    else:
        ROOT.gROOT.GetListOfCanvases().FindObject("plot").cd()
    pad = ROOT.gPad.GetPad(0)
    pad.Clear()
    pad.DivideSquare(model.indexCat().numTypes()+globs.size())
    i=0
    gConstr = {}
    if fitResult is None:
        allobs = ROOT.RooArgSet(globs); allobs.add(dataset.get())
        fitResult = ROOT.Asymptotica.makeFR(model.getParameters(allobs).selectByAttrib("Constant",False))
    for c in model.indexCat():
        i+=1
        pad.cd(i)
        ROOT.gPad.SetLeftMargin(0.2);ROOT.gPad.SetBottomMargin(0.2)
        cPdf = model.getPdf(c.first)
        for s in cPdf.servers():
            for g in globs:
                if s.dependsOn(g) and not s.dependsOn(dataset.get()):
                    gConstr[g] = s
            if s.InheritsFrom("RooRealSumPdf"):
                obs = s.getObservables(dataset).first()
                h = ROOT.TH1D("h",s.GetTitle(),obs.numBins(),obs.getBinning().array())
                h.GetXaxis().SetTitle(obs.GetTitle())
                h.GetYaxis().SetTitle("Events")
                h.SetDirectory(0);h.SetBit(ROOT.kCanDelete);ROOT.SetOwnership(h,False)
                h.SetMarkerStyle(0)
                hd = h.Clone("hdata");hd.SetMarkerStyle(20);
                ROOT.SetOwnership(hd,False)
                for j in range(0,h.GetNbinsX()):
                    obs.setBin(j)
                    h.SetBinContent(j+1, s.expectedEvents(dataset.get())*s.getVal(dataset.get())*h.GetBinWidth(j+1) )
                    h.SetBinError(j+1,ROOT.PdfWrapper.getPropagatedError2(s,fitResult,dataset.get()))
                    hd.SetBinContent(j+1, dataset.sumEntries("{}=={} && {}>={} && {}<{}".format(model.indexCat().GetName(),c.second,
                                                                          obs.GetName(),h.GetBinLowEdge(j+1),
                                                                          obs.GetName(),h.GetBinLowEdge(j+2))))
                h.Draw("hist")
                h2 = h.Clone("herr");ROOT.SetOwnership(h2,False)
                h2.SetFillStyle(3004)
                h2.SetFillColor(h2.GetLineColor())
                h2.Draw("E3same")
                hd.Draw("Esame")
                h.GetYaxis().SetRangeUser(0,max(h.GetMaximum(),hd.GetMaximum())*1.1)

    for g,c in gConstr.items():
        i+=1
        pad.cd(i)
        ROOT.gPad.SetLeftMargin(0.2);ROOT.gPad.SetBottomMargin(0.2)
        h = ROOT.TH1D("h","Aux Data",100,g.getMin(),g.getMax())
        h.GetXaxis().SetTitle(g.GetName())
        h.GetYaxis().SetTitle("Probability Density")
        h.SetDirectory(0);h.SetBit(ROOT.kCanDelete);ROOT.SetOwnership(h,False)
        for j in range(0,h.GetNbinsX()):
            c.getObservables(g).first().setVal(h.GetBinCenter(j+1))
            h.SetBinContent(j+1,c.getVal(g))
        h.Draw("C")
        l = ROOT.TLine()
        l.SetLineWidth(2)
        l.DrawLine(g.getVal(),0,g.getVal(),h.GetMaximum())
    pad.cd()
    pad.Draw()
    if oldPad: oldPad.cd()

def generate(model,fitresult,expected=False):
    snap = own(own(model.getVariables()).snapshot())
    obs = own(own(own(model.getVariables()).selectByAttrib("obs",True)).selectByAttrib("global",False))
    globs = own(own(model.getVariables()).selectByAttrib("global",True))

    if fitresult:
        own(model.getVariables()).assignValueOnly(fitresult.floatParsFinal())
        own(model.getVariables()).assignValueOnly(fitresult.constPars())

    if (expected):
        data = own(model.generate(obs, ROOT.RooFit.Extended(), ROOT.RooFit.ExpectedData(True)));
        # specific hack for this tutorial
        data_globs = own(globs.snapshot())
        for a in data_globs:
            a.setVal( own(model.getVariables()).find(a.GetName()[6:]).getVal() )
    else:
        data = own(model.generate(obs, ROOT.RooFit.Extended()))
        data_globs = own(own(model.generate(globs, 1)).get().snapshot())
    own(model.getVariables()).assignValueOnly(snap)
    return data,data_globs

def nll(model,dataset,globs):
    msglevel = ROOT.RooMsgService.instance().globalKillBelow()
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    # cannot do offsetting because will break the minNll values
    out = own(model.createNLL(dataset,ROOT.RooFit.GlobalObservables(globs),ROOT.RooFit.Offset(0)))
    ROOT.RooMsgService.instance().setGlobalKillBelow(msglevel)
    return out

def fit(model,dataset,globs):
    _nll = nll(model,dataset,globs)
    own(_nll.getVariables()).assignValueOnly(globs)
    msglevel = ROOT.RooMsgService.instance().globalKillBelow()
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    mini = ROOT.RooMinimizer(_nll)
    mini.fitter().Config().MinimizerOptions().SetPrintLevel(-1)
    mini.minimize("Minuit2")
    out = own(mini.save(ROOT.TUUID().AsString(),"Fit {} to {}".format(model.GetName(),dataset.GetName())))

    ROOT.RooMsgService.instance().setGlobalKillBelow(msglevel)
    return out
    
def asymptotic_pvalue(type,k,mu,mu_prime,sigma_mu,mu_min=-float('inf'),mu_max=float('inf')):
    return ROOT.Asymptotica.PValue(type,k,mu,mu_prime,sigma_mu,mu_min,mu_max)

def asymptotic_pvalue_qmu(qmu,mu,mu_prime,sigma_mu,mu_min=-float('inf'),mu_max=float('inf')):
    """
    Returns the asymptotic p-value of qmu test statistic 
    qmu : the one-sided-positive test statistic value 
    mu : the null hypothesis mu value
    mu_prime : the true hypothesis mu value ( = mu for null hypothesis )
    sigma_mu : stdev of mu_hat under the true hypothesis (estimate as |mu - mu_prime|/sqrt(tmu(asimov_mu_prime)))
    mu_min : min value of mu
    mu_max : max value of mu
    """
    return asymptotic_pvalue(ROOT.Asymptotica.OneSidedPositive,qmu,mu,mu_prime,sigma_mu,mu_min,mu_max)

def asymptotic_pvalue_q0(q0,mu,mu_prime,sigma_mu,mu_min=-float('inf'),mu_max=float('inf')):
    """
    q0 : the one-sided-negative test statistic value 
    mu : the null hypothesis mu value
    mu_prime : the true hypothesis mu value ( = mu for null hypothesis )
    sigma_mu : stdev of mu_hat under the true hypothesis (estimate as |mu - mu_prime|/sqrt(tmu(asimov_mu_prime))). Dependence on this parameter disappears for mu=mu_prime=0. 
    mu_min : min value of mu
    mu_max : max value of mu
    """
    return asymptotic_pvalue(ROOT.Asymptotica.OneSidedNegative,q0,mu,mu_prime,sigma_mu,mu_min,mu_max)

def asymptotic_expected_qmu(pvalue,mu,mu_prime,sigma_mu,mu_min=-float('inf'),mu_max=float('inf')):
    """
    Returns the qmu test statistic value corresponding to a given given p-value:
      pvalue : The pvalue that the returning test statistic would have under the true hypothesis (mu_prime)
      mu : The null hypothesis value of mu
      mu_prime : the true hypothesis mu value ( = mu for null hypothesis )
      sigma_mu : stdev of mu_hat under the true hypothesis (estimate as |mu - mu_prime|/sqrt(tmu(asimov_mu_prime))). Dependence on this parameter disappears for mu=mu_prime=0. 
      mu_min : min value of mu
      mu_max : max value of mu
    """
    return ROOT.Asymptotica.k(ROOT.Asymptotica.IncompatibilityFunction(ROOT.Asymptotica.OneSidedPositive,mu),
                           pvalue, mu,mu_prime,sigma_mu,mu_min,mu_max)
    
def asymptotic_expected_q0(pvalue,mu,mu_prime,sigma_mu,mu_min=-float('inf'),mu_max=float('inf')):
    """
    Returns the q0 test statistic value corresponding to a given given p-value:
      pvalue : The pvalue that the returning test statistic would have under the true hypothesis (mu_prime)
      mu : The null hypothesis value of mu
      mu_prime : the true hypothesis mu value ( = mu for null hypothesis )
      sigma_mu : stdev of mu_hat under the true hypothesis (estimate as |mu - mu_prime|/sqrt(tmu(asimov_mu_prime))). Dependence on this parameter disappears for mu=mu_prime=0. 
      mu_min : min value of mu
      mu_max : max value of mu
    """
    return ROOT.Asymptotica.k(ROOT.Asymptotica.IncompatibilityFunction(ROOT.Asymptotica.OneSidedNegative,mu),
                           pvalue, mu,mu_prime,sigma_mu,mu_min,mu_max)
    
    
def getWorkspace(day,month):
    f = ROOT.TFile("tutorialModelTemplate.root")
    w = f.Get("combined")
    
    model = w.pdf("simPdf")
    vars = ROOT.RooArgSet();
    model.treeNodeServerList(vars)
    # vars.Print("v")
    vars["mu"].setVal(1)
    data,globs = generate(model,None)
    vars["mu"].setVal(0)
    vars["sig_mass"].setVal(110)
    data.SetName("obsData");
    globs = own(own(vars.selectByAttrib("global",True)).snapshot())
    w.saveSnapshot("obsData",globs)
    w.Import(data)
    
    return w
    