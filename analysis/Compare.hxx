///////////////////////////////////////////////////////////////////////////////////
// Macro for comparing Run 2 results
///////////////////////////////////////////////////////////////////////////////////

#include "helperPlotting.h"
#include <TGraph.h>
#include <TGraphAsymmErrors.h>


void compare(const char* filenameNF, 
            const char* filenameLF,
            std::string filenameEffi,
            const std::string filenameSaving = "savingCompare.root",
            const std::string filenameSavingDirectory = "/Users/lauragans-bartl/master/MyAnalysis/EtaDalitzAnalysis/Files/"){

    ////*****FILES****/////
    TFile* fileComp = safelyOpenRootfile(filenameNF);

    //**** When input is for gammagamma
    // std::unique_ptr<TH1D> hEffiNF(fileComp->Get<TH1D>("Eta13TeV/EfficiencyEta_INT7"));
    // hEffiNF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    // hEffiNF->GetYaxis()->SetTitle("#epsilon_{NF}");
    // hEffiNF->SetName("hEffiNF");
    // setHistoStandardSettings1D(hEffiNF.get());


    // std::unique_ptr<TH1D> hEffiLF(fileComp->Get<TH1D>("Eta13TeV/EfficiencyEta_INT7B"));
    // hEffiLF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
    // hEffiLF->GetYaxis()->SetTitle("#epsilon_{LF}");
    // hEffiLF->SetName("hEffiLF");
    // setHistoStandardSettings1D(hEffiLF.get());

    //**** When input is for dalitz
    std::unique_ptr<TH1D> hEffiNF(fileComp->Get<TH1D>("TrueMesonEffiPt"));
    hEffiNF->SetName("hEffiNF");



    TFile* fileCS = safelyOpenRootfile(filenameLF);
    //**** When input is for gammagamma
    // std::unique_ptr<TGraphAsymmErrors> grEtaStat(fileCS->Get<TGraphAsymmErrors>("Eta13TeV/graphInvCrossSectionEtaComb13TeVStatErr"));
    // std::unique_ptr<TGraphAsymmErrors> grEtaSysgrEtaStat(fileCS->Get<TGraphAsymmErrors>("Eta13TeV/graphInvCrossSectionEtaComb13TeVSysErr"));
    //**** When input is for dalitz
    std::unique_ptr<TH1D> hEffiLF(fileCS->Get<TH1D>("TrueMesonEffiPt"));
    hEffiLF->SetName("hEffiLF");
    
    std::string fileEffi = filenameSavingDirectory + filenameEffi;
    TFile* fileFromEffi = safelyOpenRootfile(fileEffi);
    std::unique_ptr<TH1D> hRecEffiDalitz(fileFromEffi->Get<TH1D>("hRecEffiDalitz"));
    std::unique_ptr<TH1D> hTrueEffi(fileFromEffi->Get<TH1D>("hTrueEffi"));
    

    ////*****SAVING****/////
    std::string fileSaving = filenameSavingDirectory + filenameSaving;
    std::unique_ptr<TFile> savingFile = std::unique_ptr<TFile>(TFile::Open(fileSaving.c_str(), "RECREATE"));    //.c_str() casts std::string to char*
    hEffiNF->Write();
    hEffiLF->Write();
    hTrueEffi->Write();


    auto *hEffiNFNB = (TH1D*) hEffiLF->Clone("hEffiNFNB");
    hEffiNFNB->Reset();
    for(int newBin = 0; newBin < hEffiLF->GetXaxis()->GetNbins(); newBin++ ){
        hEffiNFNB->SetBinContent(newBin, hEffiNF->GetBinContent(newBin));
        hEffiNFNB->SetBinError(newBin, hEffiNF->GetBinError(newBin));
    }
    hEffiNFNB->SetBinContent(18, hEffiNF->GetBinContent(hEffiNF->FindBin(13)) + hEffiNF->GetBinContent(hEffiNF->FindBin(15)));
    hEffiNFNB->SetBinError(18, (hEffiNF->GetBinError(hEffiNF->FindBin(13)) * (hEffiNF->GetBinError(hEffiNF->FindBin(15)))));
    hEffiNFNB->Write();


    ////*****CALCULATING****/////
    std::vector<double> vecEffiRanges;
    for (int i = 1; i <= hEffiLF->GetNbinsX(); ++i) {
        double xE1 = hEffiLF->GetXaxis()->GetBinLowEdge(i);
        // double xE2 = hEffiLF->GetXaxis()->GetBinUpEdge(i);
        vecEffiRanges.push_back(xE1);
    }
    std::vector<double> vecPtMin = { 0., 0.2, 0.4, 0.7, 1.3, 2., 5., 18.};
    
    
    //Ratio of Run2 Effis
    auto *hEffiRatio = (TH1D*) hEffiLF->Clone("hEffiRatio");
    hEffiRatio->Divide(hEffiRatio, hEffiNFNB, 1, 1, "");       
    // hEffiRatio->SetBinContent(1, 0); // unter 0.2 ist effi unbekannt, deswegen setzen wir die einfach
    hEffiRatio->Write();
    // TGraph gr(hEffiRatio);
    // gr.Print();

    //Estimation of low field Run3 effi
    auto *hRecEffiEstiLF = (TH1D*) hTrueEffi->Clone("hRecEffiEstiLF");
    hRecEffiEstiLF->Write();
    for(int setBins = 1; setBins <= hRecEffiEstiLF->GetXaxis()->GetNbins(); setBins++){
        hRecEffiEstiLF->SetBinContent(setBins, hTrueEffi->GetBinContent(setBins) * hEffiRatio->Interpolate(hTrueEffi->GetBinCenter(setBins)));
        hRecEffiEstiLF->SetBinError(setBins, hTrueEffi->GetBinError(setBins) * hEffiRatio->Interpolate(hTrueEffi->GetBinCenter(setBins)));
        std::cout <<  setBins << ": " << hTrueEffi->GetBinContent(setBins)  << " * " << hEffiRatio->Interpolate(hTrueEffi->GetBinCenter(setBins)) << " = " << hRecEffiEstiLF->GetBinContent(setBins) << std::endl;
    }
    hRecEffiEstiLF->SetName("hRecEffiEstiLFNew");
    hRecEffiEstiLF->Write();

    std::ofstream out("debugFileComp.txt");
    for (int i = 1; i <= hEffiRatio->GetNbinsX(); ++i) {
        double x1 = hEffiRatio->GetXaxis()->GetBinLowEdge(i);
        double x2 = hEffiRatio->GetXaxis()->GetBinUpEdge(i);
        out << "x1: " << x1 << "| x2: " << x2 << "\n";
    }
    for (int i = 1; i <= hRecEffiEstiLF->GetNbinsX(); ++i) {
        double w1 = hRecEffiEstiLF->GetXaxis()->GetBinLowEdge(i);
        double q2 = hRecEffiEstiLF->GetXaxis()->GetBinUpEdge(i);
        out << "w1: " << w1 << "| q2: " << q2 << "\n";
    }

}
