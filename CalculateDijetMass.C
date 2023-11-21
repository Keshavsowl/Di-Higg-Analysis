#include <iostream>

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

void CalculateDijetMass(const char *inputFile, TH1 *histDijetMass, int color)
{
  gSystem->Load("libDelphes");

  TChain chain("Delphes");
  chain.Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  for (Long64_t entry = 0; entry < numberOfEntries; ++entry)
  {
    treeReader->ReadEntry(entry);

    int numJets = branchJet->GetEntries();

    for (int i = 0; i < numJets - 1; ++i)
    {
      Jet *jet1 = (Jet *)branchJet->At(i);

      for (int j = i + 1; j < numJets; ++j)
      {
        Jet *jet2 = (Jet *)branchJet->At(j);

        // Calculate dijet mass
        TLorentzVector lv1, lv2;
        lv1.SetPtEtaPhiM(jet1->PT, jet1->Eta, jet1->Phi, jet1->Mass);
        lv2.SetPtEtaPhiM(jet2->PT, jet2->Eta, jet2->Phi, jet2->Mass);
        TLorentzVector dijet = lv1 + lv2;
        Float_t dijetMass = dijet.M();

        histDijetMass->Fill(dijetMass);
      }
    }
  }

  // Normalize the histogram
  if (histDijetMass->Integral() > 0)
  {
    histDijetMass->Scale(1.0 / histDijetMass->Integral(), "width");
  }

  // Set the histogram color
  histDijetMass->SetLineColor(color);
}

void runDijetMassAnalysis()
{
  const char *rootFiles[] = {"HardQCD.root", "photon.root", "dihigg.root"};
  TCanvas *canvas = new TCanvas("canvas", "Dijet Mass Distribution", 800, 800);
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);

  int colors[] = {kBlue, kRed, kGreen};

  TH1 *histDijetMass[3];

  for (int i = 0; i < 3; ++i)
  {
    histDijetMass[i] = new TH1F(Form("dijetMass_%d", i), Form("Dijet Mass Distribution - %s", rootFiles[i]), 100, 0.0, 500.0); // Adjust the range as needed
    CalculateDijetMass(rootFiles[i], histDijetMass[i], colors[i]);

    histDijetMass[i]->SetMinimum(0.0);
    histDijetMass[i]->SetMaximum(0.1); // Adjust the y-axis range as needed

    if (i == 0)
      histDijetMass[i]->Draw();
    else
      histDijetMass[i]->Draw("HIST, SAME"); // Overlay histograms

    legend->AddEntry(histDijetMass[i], rootFiles[i], "l");
  }

  // Set the x-axis and y-axis labels as needed
  histDijetMass[0]->GetXaxis()->SetTitle("Dijet Mass (GeV)");
  histDijetMass[0]->GetYaxis()->SetTitle("Normalized Counts");

  legend->Draw();

  canvas->Update();
  canvas->Draw();
}

