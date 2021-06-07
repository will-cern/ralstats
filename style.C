{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.3);
  gStyle->SetHistLineWidth(1);
  gStyle->SetLineStyleString(2,"[12 12]");
  gStyle->SetEndErrorSize(0.);
  
  int font=43;
  double tsize= (font==43) ? 12 : 0.08;
  double lsize =(font==43) ? 11 : 0.07;
  gStyle->SetTextFont(font);
  gStyle->SetTextSize(tsize);
  gStyle->SetTitleOffset((font==43)? 1.5 : 0.65,"Y");
  gStyle->SetTitleFont(font,"");
  gStyle->SetTitleFontSize(tsize);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"z");
  gStyle->SetLabelSize(lsize,"x");
  gStyle->SetTitleSize(tsize,"x");
  gStyle->SetLabelSize(lsize,"y");
  gStyle->SetTitleSize(tsize,"y");
  gStyle->SetLabelSize(lsize,"z");
  gStyle->SetTitleSize(tsize,"z");
  gStyle->SetLegendFont(font);
  gStyle->SetLegendTextSize(tsize);
}