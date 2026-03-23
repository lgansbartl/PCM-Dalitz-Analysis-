///////////////////////////////////////////////////////////////////////////////////
// Macro for scanning Sgni for various pT as a function of true eta
///////////////////////////////////////////////////////////////////////////////////
#include "helperPlotting.h"
#include <cstddef>
#include <memory>


void scanSigni(const std::string& filenameSaving = "savingSigniScan.root", 
               const std::string& filenameSavingDirectory = "/Users/lauragans-bartl/master/MyAnalysis/EtaDalitzAnalysis/Files/")
    {
    
    ////*****FILES****/////
    std::vector<double> vecPtMin = { 0.1, 0.2, 0.4, 0.8, 1.2, 2., 5., 18.};
    size_t nPtSlides = vecPtMin.size();

    std::ofstream out("debugFile.txt");
    
    std::vector<std::unique_ptr<TFile>> vecfiles;
    std::vector<std::unique_ptr<TH1D>> vecHSigni;
    std::vector<std::vector<std::unique_ptr<TH1D>>> allHTrues;
    std::vector<std::vector<std::unique_ptr<TH1D>>> allHRecSignal;
    
    ////***Names***/////
    std::string filesPath = "/Users/lauragans-bartl/master/MyAnalysis/EtaDalitzAnalysis/Files/FMC/" ;
    std::vector<std::string> dirNames = {"FilesFMC100", "FilesFMC1k", "FilesFMC10k", "FilesFMC100k", "FilesFMC1M", "FilesFMC3M", "FilesFMC7M", "FilesFMC10M", "FilesFMC30M", "FilesFMC100M"};
    std::string fileName = "extrFilesMC.root";

    for(const auto& dir : dirNames){
        std::string name = filesPath + "/" + dir + "/" + fileName;
        vecfiles.push_back(std::make_unique<TFile>(name.c_str(), "READ"));
    }

    ////***Saving File***/////
    std::string fileSaving = filenameSavingDirectory + filenameSaving;
    std::unique_ptr<TFile> savingFile = std::unique_ptr<TFile>(TFile::Open(fileSaving.c_str(), "RECREATE"));
    if (!savingFile || savingFile->IsZombie()) {
        std::cerr << "Could not create output file: " << fileSaving << "\n";
        return;
    }
    
    size_t iFileIdx = 0;
    for(const auto& iFile : vecfiles){
        if (!iFile || iFile->IsZombie()) {
            std::cerr << "File invalid / not open\n";
            continue;
        }

        ////***Getting hSigni***/////
        TH1D* hSigniIn(iFile->Get<TH1D>("hSigni"));
        if (hSigniIn == nullptr) {
            std::cerr << "Histogram hSigni not found\n";
            vecHSigni.push_back(nullptr); 
            ++iFileIdx;
            continue;
        } 
        auto hSigni = std::make_unique<TH1D>(*hSigniIn);
        hSigni->SetDirectory(nullptr);                       // gut, damit es nicht an ein TFile hängt
        hSigni->SetName(Form("hSigni_%zu", iFileIdx));  
        savingFile->cd();
        hSigni->Write();
        vecHSigni.push_back(std::move(hSigni));
        ++iFileIdx;
        

        ////***Trues fuer ptSlides***/////
        std::vector<std::unique_ptr<TH1D>> histsTruePerFile;
        histsTruePerFile.reserve(vecPtMin.size() - 1);

        std::vector<std::unique_ptr<TH1D>> histsRecSignPerFile;
        histsRecSignPerFile.reserve(vecPtMin.size() - 1);

        for(size_t iSlide = 0; iSlide < vecPtMin.size()-1; iSlide++){ 
            //---Getting Trues---
            TString histNameTrue = Form("hTrue%zu", iSlide);
            TH1D* hTrueSlide = iFile->Get<TH1D>(histNameTrue);
            if (hTrueSlide==nullptr) {
                std::cerr << "histogram " << histNameTrue << " not found \n";
                continue;
            }
            auto hTrue = std::make_unique<TH1D>(*hTrueSlide);
            hTrue->SetDirectory(nullptr);
            histsTruePerFile.push_back(std::move(hTrue));
            
            //---Getting rec. Signal---
            TString histNameRec = Form("hSignal%zu", iSlide);
            TH1D* hRecSlide = iFile->Get<TH1D>(histNameRec);
            if (hRecSlide==nullptr) {
                std::cerr << "histogram " << histNameRec << " not found \n";
                continue;
            }
            auto hRecSignal = std::make_unique<TH1D>(*hRecSlide);
            hRecSignal->SetDirectory(nullptr);
            histsRecSignPerFile.push_back(std::move(hRecSignal));
        }
        allHTrues.push_back(std::move(histsTruePerFile));
        allHRecSignal.push_back(std::move(histsRecSignPerFile));
    }

 
    ////**********Pt-Slides**********/////
    double epsilon = 1.e-9;
    for(size_t iSlide = 0; iSlide < nPtSlides; iSlide++){
        auto gSigni = std::make_unique<TGraph>();
        setGraphStandardSettings1D(gSigni.get());
        gSigni->GetXaxis()->SetTitle("N_{true #eta}");
        gSigni->GetYaxis()->SetTitle("#frac{S_{rec}}{#sqrt{S_{rec}+B}}");

        auto gRecSignal = std::make_unique<TGraph>();
        setGraphStandardSettings1D(gRecSignal.get());
        gRecSignal->GetXaxis()->SetTitle("#it{N}_{true #eta}");
        gRecSignal->GetYaxis()->SetTitle("#frac{#it{N}_{#eta, rec} - #it{N}_{#eta, true}}{#it{N}_{#eta, true}}");


        out << "\n -----------iSlide-------------" <<  iSlide <<"\n";
        for(size_t iFile = 0; iFile < allHTrues.size(); iFile++ ){
            if (iFile >= vecHSigni.size() || !vecHSigni[iFile]) {
                std::cerr << "Missing vecHSigni for file " << iFile << "\n";
                continue;
            }

            ////***Calculating True Eta***/////
            if (iSlide >= allHTrues[iFile].size() || !allHTrues[iFile][iSlide]) {
                std::cerr << "Missing hTrue for file " << iFile << ", slide " << iSlide << "\n";
                continue;
            }

            std::unique_ptr<TH1D> hTrue(static_cast<TH1D*>(allHTrues[iFile][iSlide]->Clone(Form("hSame_%zu_%zu", iFile, iSlide))));
            hTrue->SetDirectory(nullptr);
            double xmin = 0.53, xmax = 0.57;
            int bmin = hTrue->GetXaxis()->FindBin(xmin);
            int bmax = hTrue->GetXaxis()->FindBin(xmax);

            out << " iFile" <<  iFile <<"\n";
            double error = 0.0;
            double nTrueEta = hTrue->IntegralAndError(bmin, bmax, error);

            double signiBin = vecHSigni[iFile]->GetXaxis()->FindBin(vecPtMin[iSlide]+epsilon);
            Double_t nSigni = vecHSigni[iFile]->GetBinContent(static_cast<Int_t>(signiBin));
            out << "# Signi: " << nSigni <<"\n \n";
            gSigni->SetPoint(static_cast<Int_t>(iFile), nTrueEta, nSigni);

            ////***Calculating Rec Eta***/////
            if (iSlide >= allHRecSignal[iFile].size() || !allHRecSignal[iFile][iSlide]) {
                std::cerr << "Missing hSignal for file " << iFile << ", slide " << iSlide << "\n";
                continue;
            }
            std::unique_ptr<TH1D> hRecSignal(static_cast<TH1D*>(allHRecSignal[iFile][iSlide]->Clone(Form("hSame_%zu_%zu", iFile, iSlide))));
            hRecSignal->SetDirectory(nullptr);

            double errorSig = 0.0;
            double nRecEta = hRecSignal->IntegralAndError(bmin, bmax, errorSig);
            // if(nTrueEta == 0)nTrueEta =1;
            if(nTrueEta == 0) continue;
            double relaError = (nRecEta - nTrueEta)/nTrueEta;
            gRecSignal->SetPoint(static_cast<Int_t>(iFile), nTrueEta, relaError);
        }
        gSigni->SetName(Form("gSigniComp_%zu", iSlide));
        gSigni->Write();

        gRecSignal->SetName(Form("gRecSignalComp_%zu", iSlide));
        gRecSignal->Write();
    }

}
