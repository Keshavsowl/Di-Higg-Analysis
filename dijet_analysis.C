#include <iostream>
#include <cmath>
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

void CalculateDeltaPhi(const char *inputFile, TH1 *histDeltaPhi, int color)
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

        // Calculate delta phi
        Float_t deltaPhi = std::abs(jet1->Phi - jet2->Phi);
        if (deltaPhi > M_PI)
          deltaPhi = 2 * M_PI - deltaPhi;

        histDeltaPhi->Fill(deltaPhi);
      }
    }
  }

  // Normalize the histogram
  if (histDeltaPhi->Integral() > 0)
  {
    histDeltaPhi->Scale(1.0 / histDeltaPhi->Integral(), "width");
  }

  // Set the histogram color
  histDeltaPhi->SetLineColor(color);
}

void CalculateDijetPt(const char *inputFile, TH1 *histDijetPt, int color)
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

        // Calculate dijet Pt
        Float_t dijetPt = (jet1->P4() + jet2->P4()).Pt();
        histDijetPt->Fill(dijetPt);
      }
    }
  }

  // Normalize the histogram
  if (histDijetPt->Integral() > 0)
  {
    histDijetPt->Scale(1.0 / histDijetPt->Integral(), "width");
  }

  // Set the histogram color
  histDijetPt->SetLineColor(color);
}

void CalculateDijetEta(const char *inputFile, TH1 *histDijetEta, int color)
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

    if (numJets >= 2)
    {
      Jet *jet1 = (Jet *)branchJet->At(0);
      Jet *jet2 = (Jet *)branchJet->At(1);

      // Calculate dijet Eta
      Float_t dijetEta = 0.5 * (jet1->Eta + jet2->Eta);
      histDijetEta->Fill(dijetEta);
    }
  }

  // Normalize the histogram
  if (histDijetEta->Integral() > 0)
  {
    histDijetEta->Scale(1.0 / histDijetEta->Integral(), "width");
  }

  // Set the histogram color
  histDijetEta->SetLineColor(color);
}

void CalculateFlavor(const char *inputFile, TH1 *histFlavor, int color)
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

        // Calculate dijet flavor
        int flavor = jet1->Flavor + jet2->Flavor;
        histFlavor->Fill(flavor);
      }
    }
  }

  // Normalize the histogram
  if (histFlavor->Integral() > 0)
  {
    histFlavor->Scale(1.0 / histFlavor->Integral(), "width");
  }

  // Set the histogram color
  histFlavor->SetLineColor(color);
}

void CalculateNCharged(const char *inputFile, TH1 *histNCharged, int color)
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

    for (int i = 0; i < numJets; ++i)
    {
      Jet *jet = (Jet *)branchJet->At(i);

      // Calculate NCharged
      int nCharged = jet->NCharged;
      histNCharged->Fill(nCharged);
    }
  }

  // Normalize the histogram
  if (histNCharged->Integral() > 0)
  {
    histNCharged->Scale(1.0 / histNCharged->Integral(), "width");
  }

  // Set the histogram color
  histNCharged->SetLineColor(color);
}

void CalculateNNeutrals(const char *inputFile, TH1 *histNNeutrals, int color)
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

    for (int i = 0; i < numJets; ++i)
    {
      Jet *jet = (Jet *)branchJet->At(i);

      // Calculate NNeutrals
      int nNeutrals = jet->NNeutrals;
      histNNeutrals->Fill(nNeutrals);
    }
  }

  // Normalize the histogram
  if (histNNeutrals->Integral() > 0)
  {
    histNNeutrals->Scale(1.0 / histNNeutrals->Integral(), "width");
  }

  // Set the histogram color
  histNNeutrals->SetLineColor(color);
}

void CalculateNeutralEnergyFraction(const char *inputFile, TH1 *histNeutralEnergyFraction, int color)
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

    for (int i = 0; i < numJets; ++i)
    {
      Jet *jet = (Jet *)branchJet->At(i);

      // Calculate NeutralEnergyFraction
      float neutralEnergyFraction = jet->NeutralEnergyFraction;
      histNeutralEnergyFraction->Fill(neutralEnergyFraction);
    }
  }

  // Normalize the histogram
  if (histNeutralEnergyFraction->Integral() > 0)
  {
    histNeutralEnergyFraction->Scale(1.0 / histNeutralEnergyFraction->Integral(), "width");
  }

  // Set the histogram color
  histNeutralEnergyFraction->SetLineColor(color);
}

void CalculateChargedEnergyFraction(const char *inputFile, TH1 *histChargedEnergyFraction, int color)
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

    for (int i = 0; i < numJets; ++i)
    {
      Jet *jet = (Jet *)branchJet->At(i);

      // Calculate ChargedEnergyFraction
      float chargedEnergyFraction = jet->ChargedEnergyFraction;
      histChargedEnergyFraction->Fill(chargedEnergyFraction);
    }
  }

  // Normalize the histogram
  if (histChargedEnergyFraction->Integral() > 0)
  {
    histChargedEnergyFraction->Scale(1.0 / histChargedEnergyFraction->Integral(), "width");
  }

  // Set the histogram color
  histChargedEnergyFraction->SetLineColor(color);
}

void runDijetAnalysisWithAdditionalHistograms()
{
  const char *rootFiles[] = {"HardQCD.root", "photon.root", "dihigg.root"};

  TCanvas *canvasDijetMass = new TCanvas("canvasDijetMass", "Dijet Mass Distribution", 800, 400);
  TCanvas *canvasDeltaPhi = new TCanvas("canvasDeltaPhi", "Delta Phi Distribution", 800, 400);
  TCanvas *canvasDijetPt = new TCanvas("canvasDijetPt", "Dijet Pt Distribution", 800, 400);
  TCanvas *canvasDijetEta = new TCanvas("canvasDijetEta", "Dijet Eta Distribution", 800, 400);
  TCanvas *canvasFlavor = new TCanvas("canvasFlavor", "Dijet Flavor Distribution", 800, 400);
  TCanvas *canvasNCharged = new TCanvas("canvasNCharged", "NCharged Distribution", 800, 400);
  TCanvas *canvasNNeutrals = new TCanvas("canvasNNeutrals", "NNeutrals Distribution", 800, 400);
  TCanvas *canvasNeutralEnergyFraction = new TCanvas("canvasNeutralEnergyFraction", "NeutralEnergyFraction Distribution", 800, 400);
  TCanvas *canvasChargedEnergyFraction = new TCanvas("canvasChargedEnergyFraction", "ChargedEnergyFraction Distribution", 800, 400);

  int colors[] = {kBlue, kRed, kGreen};

  TH1 *histDijetMass[3];
  TH1 *histDeltaPhi[3];
  TH1 *histDijetPt[3];
  TH1 *histDijetEta[3];
  TH1 *histFlavor[3];
  TH1 *histNCharged[3];
  TH1 *histNNeutrals[3];
  TH1 *histNeutralEnergyFraction[3];
  TH1 *histChargedEnergyFraction[3];

  for (int i = 0; i < 3; ++i)
  {
    // Create histograms for dijet mass, delta phi, dijet Pt, dijet Eta, flavor, NCharged, NNeutrals, NeutralEnergyFraction, and ChargedEnergyFraction
    histDijetMass[i] = new TH1F(Form("dijetMass_%d", i), Form("Dijet Mass Distribution - %s", rootFiles[i]), 100, 0.0, 500.0);
    histDeltaPhi[i] = new TH1F(Form("deltaPhi_%d", i), Form("Delta Phi Distribution - %s", rootFiles[i]), 100, 0.0, M_PI);
    histDijetPt[i] = new TH1F(Form("dijetPt_%d", i), Form("Dijet Pt Distribution - %s", rootFiles[i]), 100, 0.0, 500.0);
    histDijetEta[i] = new TH1F(Form("dijetEta_%d", i), Form("Dijet Eta Distribution - %s", rootFiles[i]), 100, -5.0, 5.0);
    histFlavor[i] = new TH1F(Form("flavor_%d", i), Form("Dijet Flavor Distribution - %s", rootFiles[i]), 7, 0.5, 6.5);
    histNCharged[i] = new TH1F(Form("nCharged_%d", i), Form("NCharged Distribution - %s", rootFiles[i]), 100, 0, 100);
    histNNeutrals[i] = new TH1F(Form("nNeutrals_%d", i), Form("NNeutrals Distribution - %s", rootFiles[i]), 100, 0, 100);
    histNeutralEnergyFraction[i] = new TH1F(Form("neutralEnergyFraction_%d", i), Form("NeutralEnergyFraction Distribution - %s", rootFiles[i]), 100, 0.0, 1.0);
    histChargedEnergyFraction[i] = new TH1F(Form("chargedEnergyFraction_%d", i), Form("ChargedEnergyFraction Distribution - %s", rootFiles[i]), 100, 0.0, 1.0);

    // Calculate dijet mass, delta phi, dijet Pt, dijet Eta, flavor, NCharged, NNeutrals, NeutralEnergyFraction, and ChargedEnergyFraction for the current file
    CalculateDijetMass(rootFiles[i], histDijetMass[i], colors[i]);
    CalculateDeltaPhi(rootFiles[i], histDeltaPhi[i], colors[i]);
    CalculateDijetPt(rootFiles[i], histDijetPt[i], colors[i]);
    CalculateDijetEta(rootFiles[i], histDijetEta[i], colors[i]);
    CalculateFlavor(rootFiles[i], histFlavor[i], colors[i]);
    CalculateNCharged(rootFiles[i], histNCharged[i], colors[i]);
    CalculateNNeutrals(rootFiles[i], histNNeutrals[i], colors[i]);
    CalculateNeutralEnergyFraction(rootFiles[i], histNeutralEnergyFraction[i], colors[i]);
    CalculateChargedEnergyFraction(rootFiles[i], histChargedEnergyFraction[i], colors[i]);

    // Set the histogram colors
    histDijetMass[i]->SetLineColor(colors[i]);
    histDeltaPhi[i]->SetLineColor(colors[i]);
    histDijetPt[i]->SetLineColor(colors[i]);
    histDijetEta[i]->SetLineColor(colors[i]);
    histFlavor[i]->SetLineColor(colors[i]);
    histNCharged[i]->SetLineColor(colors[i]);
    histNNeutrals[i]->SetLineColor(colors[i]);
    histNeutralEnergyFraction[i]->SetLineColor(colors[i]);
    histChargedEnergyFraction[i]->SetLineColor(colors[i]);

    // Normalize the histograms
    if (histDijetMass[i]->Integral() > 0)
      histDijetMass[i]->Scale(1.0 / histDijetMass[i]->Integral(), "width");
    if (histDeltaPhi[i]->Integral() > 0)
      histDeltaPhi[i]->Scale(1.0 / histDeltaPhi[i]->Integral(), "width");
    if (histDijetPt[i]->Integral() > 0)
      histDijetPt[i]->Scale(1.0 / histDijetPt[i]->Integral(), "width");
    if (histDijetEta[i]->Integral() > 0)
      histDijetEta[i]->Scale(1.0 / histDijetEta[i]->Integral(), "width");
    if (histFlavor[i]->Integral() > 0)
      histFlavor[i]->Scale(1.0 / histFlavor[i]->Integral(), "width");
    if (histNCharged[i]->Integral() > 0)
      histNCharged[i]->Scale(1.0 / histNCharged[i]->Integral(), "width");
    if (histNNeutrals[i]->Integral() > 0)
      histNNeutrals[i]->Scale(1.0 / histNNeutrals[i]->Integral(), "width");
    if (histNeutralEnergyFraction[i]->Integral() > 0)
      histNeutralEnergyFraction[i]->Scale(1.0 / histNeutralEnergyFraction[i]->Integral(), "width");
    if (histChargedEnergyFraction[i]->Integral() > 0)
      histChargedEnergyFraction[i]->Scale(1.0 / histChargedEnergyFraction[i]->Integral(), "width");

    // Set the y-axis range for each canvas
    // You can adjust this value for your specific needs
    histDijetMass[i]->SetMaximum(1.0);
    histDeltaPhi[i]->SetMaximum(1.0);
    histDijetPt[i]->SetMaximum(1.0);
    histDijetEta[i]->SetMaximum(1.0);
    histFlavor[i]->SetMaximum(1.0);
    histNCharged[i]->SetMaximum(2.0);
    histNNeutrals[i]->SetMaximum(2.0);
    histNeutralEnergyFraction[i]->SetMaximum(2.0);
    histChargedEnergyFraction[i]->SetMaximum(2.0);
  }
    // Draw dijet mass histograms on the same canvas
  canvasDijetMass->cd();
  histDijetMass[0]->Draw();
  for (int i = 1; i < 3; ++i)
    histDijetMass[i]->Draw("SAME");

  // Set axis labels for dijet mass
  histDijetMass[0]->GetXaxis()->SetTitle("Dijet Mass (GeV)");
  histDijetMass[0]->GetYaxis()->SetTitle("Normalized Counts");

  // Draw delta phi histograms on the same canvas
  canvasDeltaPhi->cd();
  histDeltaPhi[0]->Draw();
  for (int i = 1; i < 3; ++i)
    histDeltaPhi[i]->Draw("SAME");

  // Set axis labels for delta phi
  histDeltaPhi[0]->GetXaxis()->SetTitle("Delta Phi");
  histDeltaPhi[0]->GetYaxis()->SetTitle("Normalized Counts");

  // Draw dijet Pt histograms on the same canvas
  canvasDijetPt->cd();
  histDijetPt[0]->Draw();
  for (int i = 1; i < 3; ++i)
    histDijetPt[i]->Draw("SAME");

  // Set axis labels for dijet Pt
  histDijetPt[0]->GetXaxis()->SetTitle("Dijet Pt (GeV)");
  histDijetPt[0]->GetYaxis()->SetTitle("Normalized Counts");

  // Draw dijet Eta histograms on the same canvas
  canvasDijetEta->cd();
  histDijetEta[0]->Draw();
  for (int i = 1; i < 3; ++i)
    histDijetEta[i]->Draw("SAME");

  // Set axis labels for dijet Eta
  histDijetEta[0]->GetXaxis()->SetTitle("Dijet Eta");
  histDijetEta[0]->GetYaxis()->SetTitle("Normalized Counts");

  // Draw flavor histograms on the same canvas
  canvasFlavor->cd();
  histFlavor[0]->Draw();
  for (int i = 1; i < 3; ++i)
    histFlavor[i]->Draw("SAME");

  // Set axis labels for flavor
  histFlavor[0]->GetXaxis()->SetTitle("Dijet Flavor");
  histFlavor[0]->GetYaxis()->SetTitle("Normalized Counts");



  // Draw NCharged histograms on the same canvas
  canvasNCharged->cd();
  histNCharged[0]->Draw();
  for (int i = 1; i < 3; ++i)
    histNCharged[i]->Draw("SAME");

  // Set axis labels for NCharged
  histNCharged[0]->GetXaxis()->SetTitle("NCharged");
  histNCharged[0]->GetYaxis()->SetTitle("Normalized Counts");

  // Draw NNeutrals histograms on the same canvas
  canvasNNeutrals->cd();
  histNNeutrals[0]->Draw();
  for (int i = 1; i < 3; ++i)
    histNNeutrals[i]->Draw("SAME");

  // Set axis labels for NNeutrals
  histNNeutrals[0]->GetXaxis()->SetTitle("NNeutrals");
  histNNeutrals[0]->GetYaxis()->SetTitle("Normalized Counts");

  // Draw NeutralEnergyFraction histograms on the same canvas
  canvasNeutralEnergyFraction->cd();
  histNeutralEnergyFraction[0]->Draw();
  for (int i = 1; i < 3; ++i)
    histNeutralEnergyFraction[i]->Draw("SAME");

  // Set axis labels for NeutralEnergyFraction
  histNeutralEnergyFraction[0]->GetXaxis()->SetTitle("NeutralEnergyFraction");
  histNeutralEnergyFraction[0]->GetYaxis()->SetTitle("Normalized Counts");

  // Draw ChargedEnergyFraction histograms on the same canvas
  canvasChargedEnergyFraction->cd();
  histChargedEnergyFraction[0]->Draw();
  for (int i = 1; i < 3; ++i)
    histChargedEnergyFraction[i]->Draw("SAME");

  // Set axis labels for ChargedEnergyFraction
  histChargedEnergyFraction[0]->GetXaxis()->SetTitle("ChargedEnergyFraction");
  histChargedEnergyFraction[0]->GetYaxis()->SetTitle("Normalized Counts");

  // Update and display all canvases
  canvasDijetMass->Update();
  canvasDijetMass->Draw();
  canvasDeltaPhi->Update();
  canvasDeltaPhi->Draw();
  canvasDijetPt->Update();
  canvasDijetPt->Draw();
  canvasDijetEta->Update();
  canvasDijetEta->Draw();
  canvasFlavor->Update();
  canvasFlavor->Draw();
  canvasNCharged->Update();
  canvasNCharged->Draw();
  canvasNNeutrals->Update();
  canvasNNeutrals->Draw();
  canvasNeutralEnergyFraction->Update();
  canvasNeutralEnergyFraction->Draw();
  canvasChargedEnergyFraction->Update();
  canvasChargedEnergyFraction->Draw();
}


int main()
{
  runDijetAnalysisWithAdditionalHistograms();

  return 0;
}
