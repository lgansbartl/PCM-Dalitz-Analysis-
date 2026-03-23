
//******************INCLUDINGS******************//
#include "ExtractSignal.hxx"
#include "CalculateEffi.hxx"
#include "CorrectSpectrum.h"
#include "Compare.hxx"

int main(){
  
    //******************PATHS TO FILES******************//
    //***low  field 2023***//
    // const char* filename = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC23zozp.root";        //merged

    //***low  field 2025***//
    // const char* filename = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC25alV0.root";
    // const char* filename = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC25akV0.root";

    //***Nominal field 2024***//
    const char* filename = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC24.root";
    const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC24f4dR.root";
    const char* filenameMCTrue = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC24f4dT.root";

    // //***Forced Eta MC HY***//
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC26b4R.root";
    // const char* filenameMCTrue = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsLHC26b4T.root";

    //***Fast MC***//
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC100.root";
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC1k.root";
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC10k.root";
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC100k.root";
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC1M.root";
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC3M.root";
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC7M.root";
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC10M.root";
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC30M.root";
    // const char* filenameMCrec = "/Users/lauragans-bartl/Downloads/AnalysisResults/AnalysisResultsFMC100M.root";

    //***Comparing Files***//
    //* Based on Gamms Gamma
    // const char* fileRunCS = "/Users/lauragans-bartl/Downloads/AO2DAna/CombinedResultsPaperPP13TeV_2024_12_08.root";
    // const char* fileRun2 = "/Users/lauragans-bartl/Downloads/AO2DAna/InputPCMR2.root";
    //* Based on Dalitz
    const char* fileRun2NF = "/Users/lauragans-bartl/Downloads/AO2DAna/Eta_MC_GammaConv_OnlyCorrectionFactor_00010113_0dm00009f9730000dge0404000_204c6400863f02223710_0152103500000000.root";
    const char* fileRun2LF = "/Users/lauragans-bartl/Downloads/AO2DAna/Eta_MC_GammaConv_OnlyCorrectionFactor_00010113_0dm00089f9730000iih0404000_204c6400863d02263710_0152103500000000.root";



    //***File***//
    // std::string filenameSavingDirectory = std::string("/Users/lauragans-bartl/master/MyAnalysis/EtaDalitzAnalysis/Files/NominalField/FilesLHC26b4/highPt/");
    // std::string filenameSavingDirectory = std::string("/Users/lauragans-bartl/master/MyAnalysis/EtaDalitzAnalysis/Files/NominalField/FilesLHC26b4/allPt/");
    std::string filenameSavingDirectory = std::string("/Users/lauragans-bartl/master/MyAnalysis/EtaDalitzAnalysis/Files/NominalField/FilesLHC24/");


    //******************SIGNAL EXTRACTION******************//
    extractSignal(filename, "", "extrFilesData.root", filenameSavingDirectory);

    bool isMC =true;
    extractSignal(filenameMCrec, filenameMCTrue, "extrFilesMC.root", filenameSavingDirectory, isMC);

    //******************CALCULATING EFFI******************//
    std::string extrFilesMC = "extrFilesMC.root";
    std::string effiFiles = "effiFiles.root";  
    calculateEffi(extrFilesMC, effiFiles, filenameSavingDirectory);

    // //******************CORRECTING SPECTRA******************//
    //files with results from previous steps
    const char* rawYieldFilesData = "extrFilesData.root";
    const char* rawYieldFilesMC = "extrFilesMC.root";
    
    //names for files from this step
    const char* corrFilesData = "corrFilesData.root";
    const char* corrFilesMC = "corrFilesMC.root";

    correctingSpectra(rawYieldFilesData, effiFiles, corrFilesData, filenameSavingDirectory);
    correctingSpectra(rawYieldFilesMC, effiFiles, corrFilesMC, filenameSavingDirectory, isMC);
    
    // //******************COMPARRING TO RUN2 ******************//
    compare(fileRun2NF, fileRun2LF,effiFiles,  "compFilesRun2.root", filenameSavingDirectory);
}
