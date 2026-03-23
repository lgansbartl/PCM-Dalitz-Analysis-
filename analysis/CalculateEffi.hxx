///////////////////////////////////////////////////////////////////////////////////
// Calculating Raw Yield and Efficiency 
///////////////////////////////////////////////////////////////////////////////////
#include "helperPlotting.h"

void calculateEffi(const std::string& filename, 
                    const std::string filenameEffi,
                    const std::string filenameSavingDirectory = "/Users/lauragans-bartl/master/MyAnalysis/EtaDalitzAnalysis/Files/"){
    ////**********GETTING FILES**********/////
    std::string filePath = filenameSavingDirectory + filename;
    TFile* fileFromExtraction = safelyOpenRootfile(filePath); 
    std::unique_ptr<TH1D> hRawpT(fileFromExtraction->Get<TH1D>("hRawpT"));
    // hRawpT->SetDirectory(nullptr);
    std::unique_ptr<TH1D> hTrueRawpT(fileFromExtraction->Get<TH1D>("hTrueRawpT"));
    std::unique_ptr<TH1D> hGenRawpT(fileFromExtraction->Get<TH1D>("hGenRawpT"));
    std::unique_ptr<TH1D> hGenRawpTDalitz(fileFromExtraction->Get<TH1D>("hGenRawpTDalitz"));

    ////**********GDEFINITIONS**********/////
    std::string fileSaving = filenameSavingDirectory + filenameEffi;
    std::unique_ptr<TFile> effiFile = std::unique_ptr<TFile>(TFile::Open(fileSaving.c_str(), "RECREATE"));
    // hRawpT->Write();
    // hTrueRawpT->Write();
    // hGenRawpT->Write();
    // hGenRawpTDalitz->Write();

    ////**********Calculating EffICIENCY **********/////
    // Commmand to error propagation: the two hispts ar ekorrelated but not fully. When using the Divide() function, root expects the histos to be not correalted. 
    // Since they are, the errors are overestimated. The option B propagates the error binominally, wich is correct for fully correalted histos.
    //Since they are not fully correlated, using this option leads to underestimated errors
    // -> both options are not ideal!!!

    auto *hRecEffi = (TH1D*) hRawpT->Clone("hRecEffi");
    // hRecEffi->SetBinContent(0, 0.);
    // hRecEffi->SetBinError(0, 0);
    // hRecEffi->SetBinContent(1, 0.);
    // hRecEffi->SetBinError(1, 0.);
    // hRecEffi->SetBinContent(2, 0.);
    // hRecEffi->SetBinError(2, 0.);
    hRecEffi->Divide(hRecEffi, hGenRawpT.get(), 1, 1, ""); // "B"??
    setHistoStandardSettings1D(hRecEffi);
    hRecEffi->Write();

    auto *hTrueEffi = (TH1D*) hTrueRawpT->Clone("hTrueEffi");
    hTrueEffi->Divide(hTrueEffi, hGenRawpT.get(), 1, 1, ""); 
    setHistoStandardSettings1D(hTrueEffi);
    hTrueEffi->Write();

    auto *hRecEffiDalitz = (TH1D*) hRawpT->Clone("hRecEffiDalitz");
    hRecEffiDalitz->Divide(hRecEffiDalitz, hGenRawpTDalitz.get(), 1, 1, ""); 
    setHistoStandardSettings1D(hRecEffiDalitz);
    hRecEffiDalitz->Write();

    auto *hTrueEffiDalitz = (TH1D*) hTrueRawpT->Clone("hTrueEffiDalitz");
    hTrueEffiDalitz->Divide(hTrueEffiDalitz, hGenRawpTDalitz.get(), 1, 1, ""); 
    setHistoStandardSettings1D(hTrueEffiDalitz);
    hTrueEffiDalitz->Write();

}

