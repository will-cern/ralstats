
#include <cmath>
#include "Math/ProbFunc.h"
#include "Math/BrentRootFinder.h"
#include "Math/WrappedFunction.h"

#include "RooStats/RooStatsUtils.h"

class MyFitResult : public RooFitResult {
public:
    using RooFitResult::RooFitResult;
    void setCovarianceMatrix(TMatrixDSym& V) { RooFitResult::setCovarianceMatrix(V); }
    void setFinalParList(const RooArgList& l) { RooFitResult::setFinalParList(l); }
    void setInitParList(const RooArgList& l) { RooFitResult::setInitParList(l); }
    void setConstParList(const RooArgList& l) { RooFitResult::setConstParList(l); }
};

class Asymptotica {

public:

    typedef std::vector <std::pair<double, int>> IncompatFunc;

    enum PLLType {
        TwoSided = 0,
        OneSidedPositive, // for exclusions
        OneSidedNegative, // for discovery
        OneSidedAbsolute, // for exclusions by magnitude
        Uncapped, // for discovery with interest in deficits as well as excesses
        Unknown
    };

    // The incompatibility function (taking mu_hat as an input) is defined by its transitions
    // it takes values of -1, 0, or 1 ... when it 0 that means mu_hat is compatible with the hypothesis
    // Standard incompatibility functions are parameterized by mu
    // Note: the default value is taken to be 1, so an empty vector is function=1
    static IncompatFunc IncompatibilityFunction(const PLLType &type, double mu) {
        std::vector <std::pair<double, int>> out;
        if (type == TwoSided) {
            // standard PLL
        } else if (type == OneSidedPositive) {
            out.emplace_back(std::make_pair(mu, 0)); // becomes compatible @ mu_hat = mu
        } else if (type == OneSidedNegative) {
            out.emplace_back(std::make_pair(-std::numeric_limits<double>::infinity(), 0)); // compatible at -inf
            out.emplace_back(std::make_pair(mu, 1)); // becomes incompatible at mu_hat = mu
        } else if (type == OneSidedAbsolute) {
            out.emplace_back(std::make_pair(-std::numeric_limits<double>::infinity(), 0)); // compatible at -inf
            out.emplace_back(std::make_pair(-mu, 1)); // incompatible @ mu_hat = -mu
            out.emplace_back(std::make_pair(mu, 0)); // compatible again @ mu_hat = mu
        } else if (type == Uncapped) {
            out.emplace_back(std::make_pair(-std::numeric_limits<double>::infinity(), -1)); // reversed at -inf
            out.emplace_back(std::make_pair(mu, 1)); // becomes normal @ mu_hat = mu
        } else {
            throw std::runtime_error("Unknown PLL Type");
        }
        return out;
    }

    // inverse of PValue function
    static Double_t k(const IncompatFunc &compatRegions, double pValue, double poiVal, double poiPrimeVal,
                      double sigma_mu = 0,
                      double low = -std::numeric_limits<double>::infinity(),
                      double high = std::numeric_limits<double>::infinity());

    // Recommend sigma_mu = |mu - mu_prime|/sqrt(pll_mu(asimov_mu_prime))
    static Double_t PValue(const IncompatFunc &compatRegions,
                           double k,
                           double mu,
                           double mu_prime,
                           double sigma_mu = 0,
                           double mu_low = -std::numeric_limits<double>::infinity(),
                           double mu_high = std::numeric_limits<double>::infinity()
    );

    static Double_t PValue(const PLLType &pllType,
                           double k,
                           double mu,
                           double mu_prime,
                           double sigma_mu = 0,
                           double mu_low = -std::numeric_limits<double>::infinity(),
                           double mu_high = std::numeric_limits<double>::infinity()
    ) { return PValue(IncompatibilityFunction(pllType, mu), k, mu, mu_prime, sigma_mu, mu_low, mu_high); }

    static Double_t Phi_m(double mu, double mu_prime, double a, double sigma, const IncompatFunc &compatRegions);

    static int CompatFactor(const IncompatFunc &func, double mu_hat);

    static std::shared_ptr<RooFitResult> makeFR(const RooArgList& args) {

        auto fr = std::make_shared<MyFitResult>("uncorrelated");
        fr->setFinalParList(args);
        RooArgList l;
        fr->setConstParList(l);
        std::unique_ptr<RooArgList> _snap(dynamic_cast<RooArgList*>(args.snapshot()));
        fr->setInitParList(*_snap);


        TMatrixDSym cov(fr->floatParsFinal().getSize());
        int i = 0;
        for(auto& p : fr->floatParsFinal()) {
            cov(i, i) = pow(dynamic_cast<RooRealVar *>(p)->getError(), 2);
            i++;
        }
        fr->setCovarianceMatrix(cov);

        return fr;
    }

};

class PdfWrapper : public RooAbsPdf {
public:
    PdfWrapper(RooAbsPdf& f, RooAbsReal* coef=nullptr, bool expEvMode=false) : RooAbsPdf(Form("exp_%s",f.GetName())), fFunc("func","func",this,f), fCoef("coef","coef",this) {
        if (coef) fCoef.setArg(*coef);
        fExpectedEventsMode = expEvMode;
    }
    virtual ~PdfWrapper() { };
    PdfWrapper(const PdfWrapper& other, const char* name=0) : RooAbsPdf(other,name), fFunc("func",this,other.fFunc), fCoef("coef",this,other.fCoef),fExpectedEventsMode(other.fExpectedEventsMode) { }
    virtual TObject* clone(const char* newname) const override { return new PdfWrapper(*this,newname); }
    Bool_t isBinnedDistribution(const RooArgSet& obs) const override { return fFunc->isBinnedDistribution(obs); }
    std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const override {
        return fFunc->binBoundaries(obs,xlo,xhi);
    }

    double evaluate() const override {
        double out = (fExpectedEventsMode ? 1. : fFunc)*(dynamic_cast<RooAbsPdf*>(fFunc.absArg())->expectedEvents(_normSet))*(fCoef.absArg()?fCoef : 1.);
        return out;
    }

    static Double_t getPropagatedError2(RooAbsReal* func, const RooFitResult &fr, const RooArgSet &nset_in)
    {

        // Strip out parameters with zero error
        RooArgList fpf_stripped;
        RooFIter fi = fr.floatParsFinal().fwdIterator();
        RooRealVar *frv;
        while ((frv = (RooRealVar *)fi.next())) {
            if (frv->getError() > 1e-20) {
                fpf_stripped.add(*frv);
            }
        }


        // Clone self for internal use
        RooAbsReal *cloneFunc = func;//const_cast<RooAbsReal*>(dynamic_cast<const RooAbsReal*>(this)); //(RooAbsReal *)cloneTree();
        RooAbsPdf* clonePdf = dynamic_cast<RooAbsPdf*>(func);
        RooArgSet *errorParams = cloneFunc->getObservables(fpf_stripped);

        RooArgSet *nset =
                nset_in.getSize() == 0 ? cloneFunc->getParameters(*errorParams) : cloneFunc->getObservables(nset_in);

        // Make list of parameter instances of cloneFunc in order of error matrix
        RooArgList paramList;
        const RooArgList &fpf = fpf_stripped;
        vector<int> fpf_idx;
        for (Int_t i = 0; i < fpf.getSize(); i++) {
            RooAbsArg *par = errorParams->find(fpf[i].GetName());
            if (par) {
                paramList.add(*par);
                fpf_idx.push_back(i);
            }
        }

        vector<Double_t> plusVar, minusVar ;

        // Create vector of plus,minus variations for each parameter
        TMatrixDSym V(paramList.getSize()==fr.floatParsFinal().getSize()?
                      fr.covarianceMatrix():
                      fr.reducedCovarianceMatrix(paramList)) ;

        for (Int_t ivar=0 ; ivar<paramList.getSize() ; ivar++) {

            RooRealVar& rrv = (RooRealVar&)fpf[fpf_idx[ivar]] ;

            Double_t cenVal = rrv.getVal() ;
            Double_t errVal = sqrt(V(ivar,ivar)) ;

            // Make Plus variation
            ((RooRealVar*)paramList.at(ivar))->setVal(cenVal+errVal);
            plusVar.push_back(cloneFunc->getVal(nset) * (clonePdf ? clonePdf->expectedEvents(nset) : 1)) ;


            // Make Minus variation
            ((RooRealVar*)paramList.at(ivar))->setVal(cenVal-errVal) ;
            minusVar.push_back(cloneFunc->getVal(nset) * (clonePdf ? clonePdf->expectedEvents(nset) : 1)) ;

            ((RooRealVar*)paramList.at(ivar))->setVal(cenVal) ;
        }

        TMatrixDSym C(paramList.getSize()) ;
        vector<double> errVec(paramList.getSize()) ;
        for (int i=0 ; i<paramList.getSize() ; i++) {
            errVec[i] = sqrt(V(i,i)) ;
            for (int j=i ; j<paramList.getSize() ; j++) {
                C(i,j) = V(i,j)/sqrt(V(i,i)*V(j,j)) ;
                C(j,i) = C(i,j) ;
            }
        }

        // Make vector of variations
        TVectorD F(plusVar.size()) ;
        for (unsigned int j=0 ; j<plusVar.size() ; j++) {
            F[j] = (plusVar[j]-minusVar[j])/2 ;
        }

        // Calculate error in linear approximation from variations and correlation coefficient
        Double_t sum = F*(C*F) ;

        //delete cloneFunc ;
        delete errorParams ;
        delete nset ;


        return sqrt(sum) ;
    }


private:
    RooRealProxy fFunc;
    RooRealProxy fCoef;
    bool fExpectedEventsMode=false;
};




Double_t Asymptotica::k(const IncompatFunc& compatRegions, double pValue, double poiVal, double poiPrimeVal, double sigma, double low, double high) {

    // determine the pll value corresponding to nSigma expected - i.e. where the altPValue equals e.g. 50% for nSigma=0,
    //find the solution (wrt x) of: FitManager::altPValue(x, var(poi), alt_val, _sigma_mu, _compatibilityFunction) - targetPValue = 0
    double targetTailIntegral = pValue; //ROOT::Math::normal_cdf(*nSigma);

    // check how much of the alt distribution density is in the delta function @ 0
    // if more than 1 - target is in there, if so then pll must be 0
    double prob_in_delta = Phi_m(poiVal, poiPrimeVal, std::numeric_limits<double>::infinity(), sigma, compatRegions);
    // also get a contribution to the delta function for mu_hat < mu_L IF mu==mu_L
    if (poiVal == low) {
        // since mu_hat is gaussian distributed about mu_prime with std-dev = sigma
        // the integral is Phi( mu_L - mu_prime / (sigma) )
        double mu_L = low;
        prob_in_delta += ROOT::Math::normal_cdf((mu_L - poiPrimeVal) / sigma);
    }

    if (prob_in_delta > 1-targetTailIntegral) {
        return 0;
    }

    struct TailIntegralFunction {
        TailIntegralFunction(double _poiVal, double _alt_val, double _sigma_mu, double _low, double _high,
                             IncompatFunc _compatibilityFunction, double _target) :
                poiVal(_poiVal), alt_val(_alt_val), sigma_mu(_sigma_mu), low(_low),high(_high), target(_target),
                cFunc(_compatibilityFunction) {}

        double operator()(double x) const {
            return PValue(cFunc, x, poiVal, alt_val, sigma_mu,low,high) - target;
        }

        double poiVal, alt_val, sigma_mu, low,high,target;
        IncompatFunc cFunc;
    };

    TailIntegralFunction f(poiVal, poiPrimeVal, sigma, low, high, compatRegions, targetTailIntegral);
    ROOT::Math::BrentRootFinder brf;
    ROOT::Math::WrappedFunction<TailIntegralFunction> wf(f);

    auto tmpLvl = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;
    double _pll = 500.;
    double currVal(1.);
    int tryCount(0);
    double _prev_pll = _pll;
    do {
        currVal = wf(_pll);
        if (currVal > 1e-4) _pll = 2. * (_pll + 1.); // goto bigger pll scale
        else if (currVal < -1e-4) _pll /= 2.; // goto smaller pll scale
        //std::cout << "pll = " << _pll << " currVal = " << currVal << std::endl;
        brf.SetFunction(wf, 0, _pll);
        if (brf.Solve()) {
            _prev_pll = _pll;
            _pll = brf.Root();
        }
        //std::cout << " -- " << brf.Root() << " " << FitManager::altPValue(_pll, mu, alt_val, sigma, pllModifier()) << " >> " << wf(_pll) << std::endl;
        tryCount++;
        if (tryCount > 20) {
            gErrorIgnoreLevel = tmpLvl;
            Warning("Asymptotica::k", "Reached limit pValue=%g pll=%g", pValue, _pll);
            break;
        }
    } while (std::abs(wf(_pll)) > 1e-4 && std::abs(wf(_pll)) < std::abs(wf(_prev_pll)) * 0.99);
    gErrorIgnoreLevel = tmpLvl;
    return _pll;

}

Double_t Asymptotica::PValue(const IncompatFunc& compatRegions, double k, double poiVal, double poi_primeVal, double sigma, double lowBound, double upBound) {
    //uncapped test statistic is equal to onesidednegative when k is positive, and equal to 1.0 - difference between twosided and onesidednegative when k is negative ...
    if(compatRegions==IncompatibilityFunction(Uncapped,poiVal)) {
        //if(k==0) return 0.5;
        if(k>0) return PValue(OneSidedNegative,k,poiVal,poi_primeVal,sigma,lowBound,upBound);
        return 1. - (PValue(TwoSided,-k,poiVal,poi_primeVal,sigma,lowBound,upBound) - PValue(OneSidedNegative,-k,poiVal,poi_primeVal,sigma,lowBound,upBound));
    }


    //if(k<0) return 1.;
    if(k<=0) {
        if(compatRegions==IncompatibilityFunction(OneSidedNegative,poiVal) && std::abs(poiVal-poi_primeVal)<1e-9) return 0.5; //when doing discovery (one-sided negative) use a 0.5 pValue
        return 1.; //case to catch the delta function that ends up at exactly 0 for the one-sided tests
    }

    if(sigma==0) {
        if (lowBound!=-std::numeric_limits<double>::infinity() || upBound != std::numeric_limits<double>::infinity()) {
            return -1;
        } else if(std::abs(poiVal-poi_primeVal)>1e-12) {
            return -1;
        }
    }

    //get the poi value that defines the test statistic, and the poi_prime hypothesis we are testing
    //when setting limits, these are often the same value

    Double_t Lambda_y = 0;
    if(std::abs(poiVal-poi_primeVal)>1e-12) Lambda_y = (poiVal-poi_primeVal)/sigma;

    Double_t k_low = (lowBound == -std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::infinity() : pow((poiVal - lowBound)/sigma,2);
    Double_t k_high = (upBound == std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::infinity() : pow((upBound - poiVal)/sigma,2);

    double out = Phi_m(poiVal,poi_primeVal,std::numeric_limits<double>::infinity(),sigma,compatRegions) - 1;

    if (out < -1) {
        // compatibility function is unsupported, return negative
        return -2;
    }

    //go through the 4 'regions' ... only two of which will apply
    if( k <= k_high ) {
        out += ROOT::Math::gaussian_cdf(sqrt(k)+Lambda_y) - Phi_m(poiVal,poi_primeVal,Lambda_y + sqrt(k),sigma,compatRegions);
    } else {
        double Lambda_high = (poiVal - upBound)*(poiVal + upBound - 2.*poi_primeVal)/(sigma*sigma);
        double sigma_high = 2.*(upBound-poiVal)/sigma;
        out +=  ROOT::Math::gaussian_cdf((k-Lambda_high)/sigma_high) - Phi_m(poiVal,poi_primeVal,(k - Lambda_high)/sigma_high,sigma,compatRegions);
    }

    if( k <= k_low ) {
        out += ROOT::Math::gaussian_cdf(sqrt(k)-Lambda_y) + Phi_m(poiVal,poi_primeVal,Lambda_y - sqrt(k),sigma,compatRegions);
    } else {
        double Lambda_low = (poiVal - lowBound)*(poiVal + lowBound - 2.*poi_primeVal)/(sigma*sigma);
        double sigma_low = 2.*(poiVal - lowBound)/sigma;
        out +=  ROOT::Math::gaussian_cdf((k-Lambda_low)/sigma_low) + Phi_m(poiVal,poi_primeVal,(Lambda_low - k)/sigma_low,sigma,compatRegions);
        /*out +=  ROOT::Math::gaussian_cdf((k-Lambda_low)/sigma_low) +
             2*Phi_m(poiVal,poi_primeVal,(Lambda_low - k_low)==0 ? 0 : ((Lambda_low - k_low)/sigma_low),sigma,compatRegions)
             - Phi_m(poiVal,poi_primeVal,(Lambda_low - k)/sigma_low,sigma,compatFunc);
*/

        // handle case where poiVal = lowBound (e.g. testing mu=0 when lower bound is mu=0).
        // sigma_low will be 0 and gaussian_cdf will end up being 1, but we need it to converge instead
        // to 0.5 so that pValue(k=0) converges to 1 rather than 0.5.
        // handle this by 'adding' back in the lower bound
        // TODO: Think more about this?
        /*if (sigma_low == 0) {
            out -= 0.5;
        }*/


    }

    return 1. - out;
}



Double_t Asymptotica::Phi_m(double mu, double mu_prime, double a, double sigma, const IncompatFunc& compatRegions ) {

    if (sigma==0) sigma=1e-100;// avoid nans if sigma is 0

    // want to evaluate gaussian integral in regions where IncompatFunc = 0

    double out = 0;
    double lowEdge = std::numeric_limits<double>::quiet_NaN();
    for(auto& transition : compatRegions) {
        if (transition.first >= (a*sigma + mu_prime)) break;
        if(transition.second==0 && std::isnan(lowEdge)) lowEdge = transition.first;
        else if(!std::isnan(lowEdge)) {
            out += ROOT::Math::gaussian_cdf((transition.first-mu_prime)/sigma) - ROOT::Math::gaussian_cdf((lowEdge-mu_prime)/sigma);
            lowEdge = std::numeric_limits<double>::quiet_NaN();
        }
    }
    if (!std::isnan(lowEdge)) { // also catches case where last transition is before a
        out += ROOT::Math::gaussian_cdf(a) - ROOT::Math::gaussian_cdf((lowEdge-mu_prime)/sigma);
    }

    return out;

}

int Asymptotica::CompatFactor(const IncompatFunc& func, double mu_hat) {
    if (std::isnan(mu_hat)) return 1; // nan is never compatible
    int out = 1;
    for(auto& transition : func) {
        if (transition.first > mu_hat) break;
        out = transition.second;
    }
    return out;
}
