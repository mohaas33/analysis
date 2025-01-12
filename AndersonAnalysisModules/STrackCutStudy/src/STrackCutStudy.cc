// ----------------------------------------------------------------------------
// 'STrackCutStudy.cc'
// Derek Anderson
// 12.15.2022
//
// Reads in the 'ntp_track' Ntuple
// generated by the SVtxEvaluator
// class and studies the impact
// of cutting on various quantities.
// ----------------------------------------------------------------------------

#define STRACKCUTSTUDY_CC

// header files
#include "STrackCutStudy.h"
#include "STrackCutStudy.io.h"
#include "STrackCutStudy.ana.h"
#include "STrackCutStudy.hist.h"
#include "STrackCutStudy.plot.h"

using namespace std;



// ctor/dtor ------------------------------------------------------------------

STrackCutStudy::STrackCutStudy() {

  // clear member variables
  sTxtEO.clear();
  sTxtPU.clear();
  nTxtEO          = 0;
  nTxtPU          = 0;
  inBatchMode     = false;
  makePlots       = false;
  doPileup        = false;
  doIntNorm       = false;
  doBeforeCuts    = false;
  doAvgClustCalc  = false;
  normalPtFracMin = 0.;
  normalPtFracMax = 9999.;
  doPrimaryCut    = false;
  doMVtxCut       = false;
  doVzCut         = false;
  doDcaXyCut      = false;
  doDcaZCut       = false;
  doQualityCut    = false;

  // set whether or not type is truth
  isTruth[TYPE::TRACK]         = false;
  isTruth[TYPE::TRUTH]         = true;
  isTruth[TYPE::WEIRD_ALL]     = false;
  isTruth[TYPE::WEIRD_SI]      = false;
  isTruth[TYPE::WEIRD_TPC]     = false;
  isTruth[TYPE::NORMAL]        = false;
  isTruth[TYPE::PILEUP]        = false;
  isTruth[TYPE::PRIMARY]       = false;
  isTruth[TYPE::NONPRIM]       = false;
  isTruth[TYPE::TRK_CUT]       = false;
  isTruth[TYPE::TRU_CUT]       = true;
  isTruth[TYPE::WEIRD_CUT]     = false;
  isTruth[TYPE::WEIRD_SI_CUT]  = false;
  isTruth[TYPE::WEIRD_TPC_CUT] = false;
  isTruth[TYPE::NORM_CUT]      = false;
  isTruth[TYPE::PILE_CUT]      = false;
  isTruth[TYPE::PRIM_CUT]      = false;
  isTruth[TYPE::NONPRIM_CUT]   = false;

  // set whether or not type has pileup
  isPileup[TYPE::TRACK]         = false;
  isPileup[TYPE::TRUTH]         = false;
  isPileup[TYPE::WEIRD_ALL]     = false;
  isPileup[TYPE::WEIRD_SI]      = false;
  isPileup[TYPE::WEIRD_TPC]     = false;
  isPileup[TYPE::NORMAL]        = false;
  isPileup[TYPE::PILEUP]        = true;
  isPileup[TYPE::PRIMARY]       = true;
  isPileup[TYPE::NONPRIM]       = true;
  isPileup[TYPE::TRK_CUT]       = false;
  isPileup[TYPE::TRU_CUT]       = false;
  isPileup[TYPE::WEIRD_CUT]     = false;
  isPileup[TYPE::WEIRD_SI_CUT]  = false;
  isPileup[TYPE::WEIRD_TPC_CUT] = false;
  isPileup[TYPE::NORM_CUT]      = false;
  isPileup[TYPE::PILE_CUT]      = true;
  isPileup[TYPE::PRIM_CUT]      = true;
  isPileup[TYPE::NONPRIM_CUT]   = true;

  // set whether or not type is before cuts
  isBeforeCuts[TYPE::TRACK]         = true;
  isBeforeCuts[TYPE::TRUTH]         = true;
  isBeforeCuts[TYPE::WEIRD_ALL]     = true;
  isBeforeCuts[TYPE::WEIRD_SI]      = true;
  isBeforeCuts[TYPE::WEIRD_TPC]     = true;
  isBeforeCuts[TYPE::NORMAL]        = true;
  isBeforeCuts[TYPE::PILEUP]        = true;
  isBeforeCuts[TYPE::PRIMARY]       = true;
  isBeforeCuts[TYPE::NONPRIM]       = true;
  isBeforeCuts[TYPE::TRK_CUT]       = false;
  isBeforeCuts[TYPE::TRU_CUT]       = false;
  isBeforeCuts[TYPE::WEIRD_CUT]     = false;
  isBeforeCuts[TYPE::WEIRD_SI_CUT]  = false;
  isBeforeCuts[TYPE::WEIRD_TPC_CUT] = false;
  isBeforeCuts[TYPE::NORM_CUT]      = false;
  isBeforeCuts[TYPE::PILE_CUT]      = false;
  isBeforeCuts[TYPE::PRIM_CUT]      = false;
  isBeforeCuts[TYPE::NONPRIM_CUT]   = false;

  // set whether or not track variable has a truth value
  trkVarHasTruVal[TRKVAR::VX]       = true;
  trkVarHasTruVal[TRKVAR::VY]       = true;
  trkVarHasTruVal[TRKVAR::VZ]       = true;
  trkVarHasTruVal[TRKVAR::NMMS]     = true;
  trkVarHasTruVal[TRKVAR::NMAP]     = true;
  trkVarHasTruVal[TRKVAR::NINT]     = true;
  trkVarHasTruVal[TRKVAR::NTPC]     = true;
  trkVarHasTruVal[TRKVAR::QUAL]     = false;
  trkVarHasTruVal[TRKVAR::DCAXY]    = false;
  trkVarHasTruVal[TRKVAR::DCAZ]     = false;
  trkVarHasTruVal[TRKVAR::DELDCAXY] = false;
  trkVarHasTruVal[TRKVAR::DELDCAZ]  = false;
  trkVarHasTruVal[TRKVAR::NCLUST]   = true;
  trkVarHasTruVal[TRKVAR::AVGCLUST] = true;

  // set whether or not physics variable has a truth value
  physVarHasTruVal[PHYSVAR::PHI]    = true;
  physVarHasTruVal[PHYSVAR::ETA]    = true;
  physVarHasTruVal[PHYSVAR::PT]     = true;
  physVarHasTruVal[PHYSVAR::DELPHI] = false;
  physVarHasTruVal[PHYSVAR::DELETA] = false;
  physVarHasTruVal[PHYSVAR::DELPT]  = false;

  // set type colors
  fTypeCol[TYPE::TRACK]         = 923;
  fTypeCol[TYPE::TRUTH]         = 899;
  fTypeCol[TYPE::WEIRD_ALL]     = 879;
  fTypeCol[TYPE::WEIRD_SI]      = 809;
  fTypeCol[TYPE::WEIRD_TPC]     = 849;
  fTypeCol[TYPE::NORMAL]        = 889;
  fTypeCol[TYPE::PILEUP]        = 923;
  fTypeCol[TYPE::PRIMARY]       = 859;
  fTypeCol[TYPE::NONPRIM]       = 799;
  fTypeCol[TYPE::TRK_CUT]       = 923;
  fTypeCol[TYPE::TRU_CUT]       = 899;
  fTypeCol[TYPE::WEIRD_CUT]     = 879;
  fTypeCol[TYPE::WEIRD_SI_CUT]  = 809;
  fTypeCol[TYPE::WEIRD_TPC_CUT] = 849;
  fTypeCol[TYPE::NORM_CUT]      = 889;
  fTypeCol[TYPE::PILE_CUT]      = 923;
  fTypeCol[TYPE::PRIM_CUT]      = 859;
  fTypeCol[TYPE::NONPRIM_CUT]   = 799;

  // set type markers
  fTypeMar[TYPE::TRACK]         = 20;
  fTypeMar[TYPE::TRUTH]         = 24;
  fTypeMar[TYPE::WEIRD_ALL]     = 26;
  fTypeMar[TYPE::WEIRD_SI]      = 5;
  fTypeMar[TYPE::WEIRD_TPC]     = 2;
  fTypeMar[TYPE::NORMAL]        = 32;
  fTypeMar[TYPE::PILEUP]        = 20;
  fTypeMar[TYPE::PRIMARY]       = 26;
  fTypeMar[TYPE::NONPRIM]       = 32;
  fTypeMar[TYPE::TRK_CUT]       = 20;
  fTypeMar[TYPE::TRU_CUT]       = 24;
  fTypeMar[TYPE::WEIRD_CUT]     = 26;
  fTypeMar[TYPE::WEIRD_SI_CUT]  = 5;
  fTypeMar[TYPE::WEIRD_TPC_CUT] = 2;
  fTypeMar[TYPE::NORM_CUT]      = 32;
  fTypeMar[TYPE::PILE_CUT]      = 20;
  fTypeMar[TYPE::PRIM_CUT]      = 26;
  fTypeMar[TYPE::NONPRIM_CUT]   = 32;

  // set type names
  sTrkNames[TYPE::TRACK]         = "AllTrack";
  sTrkNames[TYPE::TRUTH]         = "AllTruth";
  sTrkNames[TYPE::WEIRD_ALL]     = "AllWeird";
  sTrkNames[TYPE::WEIRD_SI]      = "AllSiWeird";
  sTrkNames[TYPE::WEIRD_TPC]     = "AllTpcWeird";
  sTrkNames[TYPE::NORMAL]        = "AllNormal";
  sTrkNames[TYPE::PILEUP]        = "AllPileup";
  sTrkNames[TYPE::PRIMARY]       = "AllPrimePileup";
  sTrkNames[TYPE::NONPRIM]       = "AllNonPrimePileup";
  sTrkNames[TYPE::TRK_CUT]       = "CutTrack";
  sTrkNames[TYPE::TRU_CUT]       = "CutTruth";
  sTrkNames[TYPE::WEIRD_CUT]     = "CutWeird";
  sTrkNames[TYPE::WEIRD_SI_CUT]  = "CutSiWeird";
  sTrkNames[TYPE::WEIRD_TPC_CUT] = "CutTpcWeird";
  sTrkNames[TYPE::NORM_CUT]      = "CutNormal";
  sTrkNames[TYPE::PILE_CUT]      = "CutPileup";
  sTrkNames[TYPE::PRIM_CUT]      = "CutPrimePileup";
  sTrkNames[TYPE::NONPRIM_CUT]   = "CutNonPrimePileup"; 

  // set type plot labels
  sTrkLabels[TYPE::TRACK]         = "Tracks (before cuts)";
  sTrkLabels[TYPE::TRUTH]         = "Truth tracks (before cuts)";
  sTrkLabels[TYPE::WEIRD_ALL]     = "Weird tracks (before cuts)";
  sTrkLabels[TYPE::WEIRD_SI]      = "Weird tracks (Si seed, before cuts)";
  sTrkLabels[TYPE::WEIRD_TPC]     = "Weird tracks (TPC seed, before cuts)";
  sTrkLabels[TYPE::NORMAL]        = "Normal tracks (before cuts)";
  sTrkLabels[TYPE::PILEUP]        = "Including pileup tracks (all, before cuts)";
  sTrkLabels[TYPE::PRIMARY]       = "Including pileup tracks (only primary, before cuts)";
  sTrkLabels[TYPE::NONPRIM]       = "Including pileup gracks (non-primary, before cuts)";
  sTrkLabels[TYPE::TRK_CUT]       = "Tracks (after cuts)";
  sTrkLabels[TYPE::TRU_CUT]       = "Truth tracks (after cuts)";
  sTrkLabels[TYPE::WEIRD_CUT]     = "Weird tracks (after cuts)";
  sTrkLabels[TYPE::WEIRD_SI_CUT]  = "Weird tracks (Si seed, after cuts)";
  sTrkLabels[TYPE::WEIRD_TPC_CUT] = "Weird tracks (TPC seed, after cuts)";
  sTrkLabels[TYPE::NORM_CUT]      = "Normal tracks (after cuts)";
  sTrkLabels[TYPE::PILE_CUT]      = "Including pileup tracks (all, after cuts)";
  sTrkLabels[TYPE::PRIM_CUT]      = "Including pileup tracks (only primary, after cuts)";
  sTrkLabels[TYPE::NONPRIM_CUT]   = "Including pileup gracks (non-primary, after cuts)";

  // set track variable names
  sTrkVars[TRKVAR::VX]       = "Vx";
  sTrkVars[TRKVAR::VY]       = "Vy";
  sTrkVars[TRKVAR::VZ]       = "Vz";
  sTrkVars[TRKVAR::NMMS]     = "NMms";
  sTrkVars[TRKVAR::NMAP]     = "NMap";
  sTrkVars[TRKVAR::NINT]     = "NInt";
  sTrkVars[TRKVAR::NTPC]     = "NTpc";
  sTrkVars[TRKVAR::QUAL]     = "Qual";
  sTrkVars[TRKVAR::DCAXY]    = "DcaXY";
  sTrkVars[TRKVAR::DCAZ]     = "DcaZ";
  sTrkVars[TRKVAR::DELDCAXY] = "DeltaDcaXY";
  sTrkVars[TRKVAR::DELDCAZ]  = "DeltaDcaZ";
  sTrkVars[TRKVAR::NCLUST]   = "NClust";
  sTrkVars[TRKVAR::AVGCLUST] = "AvgClustSize";

  // set physics variable names
  sPhysVars[PHYSVAR::PHI]    = "Phi";
  sPhysVars[PHYSVAR::ETA]    = "Eta";
  sPhysVars[PHYSVAR::PT]     = "Pt";
  sPhysVars[PHYSVAR::DELPHI] = "DeltaPhi";
  sPhysVars[PHYSVAR::DELETA] = "DeltaEta";
  sPhysVars[PHYSVAR::DELPT]  = "DeltaPt";
  cout << "\n  Beginning track cut study."  << endl;

}  // end ctor




STrackCutStudy::~STrackCutStudy() {

  const bool doTuplesExist = (ntTrkEO || ntTrkPU);
  if (!doTuplesExist) {
    return;
  } else {
    delete ntTrkEO -> GetCurrentFile();
    delete ntTrkPU -> GetCurrentFile();
  }

}  // end dtor



//  main methods --------------------------------------------------------------

void STrackCutStudy::Init() {

  // announce method
  cout << "    Initializing:" << endl;
  InitFiles();
  InitTuples();
  InitHists();
  MakeCutText();
  return;

}  // end Init()



void STrackCutStudy::Analyze() {

  // check for tree and announce method
  const bool doTuplesExist = (ntTrkEO && ntTrkPU);
  if (!doTuplesExist) {
    cerr << "PANIC: no input tuples!\n"
         << "        ntTrkEO = " << ntTrkEO << ", ntTrkPU = " << ntTrkPU
         << endl;
    assert(doTuplesExist);
  }
  cout << "    Analyzing:" <<endl;

  // prepare for embed-only entry loop
  Long64_t nEntriesEO = ntTrkEO -> GetEntries();
  cout << "      Beginning embed-only entry loop: " << nEntriesEO << " entries to process..." << endl;

  // arrays for filling histograms
  Double_t recoTrkVars[NTrkVar];
  Double_t trueTrkVars[NTrkVar];
  Double_t recoPhysVars[NPhysVar];
  Double_t truePhysVars[NPhysVar];

  // loop over embed-only tuple entries
  Long64_t nBytesEO(0);
  for (Long64_t iEntry = 0; iEntry < nEntriesEO; iEntry++) {

    // grab entry
    const Long64_t bytesEO = ntTrkEO -> GetEntry(iEntry);
    if (bytesEO < 0.) {
      cerr << "WARNING: something wrong with embed-only entry #" << iEntry << "! Aborting loop!" << endl;
      break;
    }
    nBytesEO += bytesEO;

    // announce progress
    const Long64_t iProg = iEntry + 1;
    if (inBatchMode) {
      cout << "        Processing embed-only entry " << iProg << "/" << nEntriesEO << "..." << endl;
    } else {
      if (iProg == nEntriesEO) {
        cout << "        Processing embed-only entry " << iProg << "/" << nEntriesEO << "..." << endl;
      } else {
        cout << "        Processing embed-only entry " << iProg << "/" << nEntriesEO << "...\r" << flush;
      }
    }

    // perform calculations
    const Double_t umDcaXY    = dca3dxy * 10000;
    const Double_t umDcaZ     = dca3dz * 10000;
    const Double_t deltaDcaXY = abs(dca3dxysigma / dca3dxy);
    const Double_t deltaDcaZ  = abs(dca3dzsigma / dca3dz);
    const Double_t deltaEta   = abs(deltaeta / eta);
    const Double_t deltaPhi   = abs(deltaphi / phi);
    const Double_t deltaPt    = abs(deltapt / pt);
    const Double_t ptFrac     = pt /gpt;

    // set reco track variables
    recoTrkVars[TRKVAR::VX]       = vx;
    recoTrkVars[TRKVAR::VY]       = vy;
    recoTrkVars[TRKVAR::VZ]       = vz;
    recoTrkVars[TRKVAR::NMMS]     = (Double_t) nlmms;
    recoTrkVars[TRKVAR::NMAP]     = (Double_t) nlmaps;
    recoTrkVars[TRKVAR::NINT]     = (Double_t) nlintt;
    recoTrkVars[TRKVAR::NTPC]     = (Double_t) ntpc;
    recoTrkVars[TRKVAR::QUAL]     = quality;
    recoTrkVars[TRKVAR::DCAXY]    = umDcaXY;
    recoTrkVars[TRKVAR::DCAZ]     = umDcaZ;
    recoTrkVars[TRKVAR::DELDCAXY] = deltaDcaXY;
    recoTrkVars[TRKVAR::DELDCAZ]  = deltaDcaZ;
    recoTrkVars[TRKVAR::NCLUST]   = 0.;
    recoTrkVars[TRKVAR::AVGCLUST] = 0.;

    // set true track variables
    trueTrkVars[TRKVAR::VX]       = gvx;
    trueTrkVars[TRKVAR::VY]       = gvy;
    trueTrkVars[TRKVAR::VZ]       = gvz;
    trueTrkVars[TRKVAR::NMMS]     = (Double_t) gnlmms;
    trueTrkVars[TRKVAR::NMAP]     = (Double_t) gnlmaps;
    trueTrkVars[TRKVAR::NINT]     = (Double_t) gnlintt;
    trueTrkVars[TRKVAR::NTPC]     = (Double_t) gntpc;
    trueTrkVars[TRKVAR::QUAL]     = quality;
    trueTrkVars[TRKVAR::DCAXY]    = umDcaXY;
    trueTrkVars[TRKVAR::DCAZ]     = umDcaZ;
    trueTrkVars[TRKVAR::DELDCAXY] = deltaDcaXY;
    trueTrkVars[TRKVAR::DELDCAZ]  = deltaDcaZ;
    trueTrkVars[TRKVAR::NCLUST]   = 0.;
    trueTrkVars[TRKVAR::AVGCLUST] = 0.;

    // set reco phys variables
    recoPhysVars[PHYSVAR::PHI]    = phi;
    recoPhysVars[PHYSVAR::ETA]    = eta;
    recoPhysVars[PHYSVAR::PT]     = pt;
    recoPhysVars[PHYSVAR::DELPHI] = deltaPhi;
    recoPhysVars[PHYSVAR::DELETA] = deltaEta;
    recoPhysVars[PHYSVAR::DELPT]  = deltaPt;

    // set true phys variables
    truePhysVars[PHYSVAR::PHI]    = gphi;
    truePhysVars[PHYSVAR::ETA]    = geta;
    truePhysVars[PHYSVAR::PT]     = gpt;
    truePhysVars[PHYSVAR::DELPHI] = deltaPhi;
    truePhysVars[PHYSVAR::DELETA] = deltaEta;
    truePhysVars[PHYSVAR::DELPT]  = deltaPt;

    // check for weird tracks
    const Bool_t hasSiSeed    = (nmaps == 3);
    const Bool_t hasTpcSeed   = (nmaps == 0);
    const Bool_t isPrimary    = (gprimary == 1);
    const Bool_t isWeirdTrack = ((ptFrac < normalPtFracMin) || (ptFrac > normalPtFracMax));

    // fill all track histograms
    if (doBeforeCuts) {
      FillTrackHistograms(TYPE::TRACK, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
      FillTruthHistograms(TYPE::TRUTH, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
 
      // fill all embed_only weird histograms
      if (isWeirdTrack) {
        FillTrackHistograms(TYPE::WEIRD_ALL, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
        if (hasSiSeed)  FillTrackHistograms(TYPE::WEIRD_SI,  recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
        if (hasTpcSeed) FillTrackHistograms(TYPE::WEIRD_TPC, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
      } else {
        FillTrackHistograms(TYPE::NORMAL, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
      }
    }

    // apply cuts
    const Bool_t isGoodTrk = ApplyCuts(isPrimary, (UInt_t) nlmaps, (UInt_t) ntpc, vz, umDcaXY, umDcaZ, quality);
    if (!isGoodTrk) continue;

    // fill cut track histograms
    FillTrackHistograms(TYPE::TRK_CUT, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
    FillTruthHistograms(TYPE::TRU_CUT, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
 
    // fill cut embed_only weird histograms
    if (isWeirdTrack) {
      FillTrackHistograms(TYPE::WEIRD_CUT, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
      if (hasSiSeed)  FillTrackHistograms(TYPE::WEIRD_SI_CUT,  recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
      if (hasTpcSeed) FillTrackHistograms(TYPE::WEIRD_TPC_CUT, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
    } else {
      FillTrackHistograms(TYPE::NORM_CUT, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
    }
  }  // end embed-only entry loop
  cout << "      Finished embed-only entry loop." << endl;

  // prepare for with-pileup entry loop
  if (doPileup) {
    Long64_t nEntriesPU = ntTrkPU -> GetEntries();
    cout << "      Beginning with-pileup entry loop: " << nEntriesPU << " entries to process..." << endl;

    // loop over with-pileup tuple entries
    Long64_t nBytesPU(0);
    for (Long64_t iEntry = 0; iEntry < nEntriesPU; iEntry++) {

      // grab entry
      const Long64_t bytesPU = ntTrkPU -> GetEntry(iEntry);
      if (bytesPU < 0.) {
        cerr << "WARNING: something wrong with with-pileup entry #" << iEntry << "! Aborting loop!" << endl;
        break;
      }
      nBytesPU += bytesPU;

      // announce progress
      const Long64_t iProg = iEntry + 1;
      if (inBatchMode) {
        cout << "        Processing with-pileup entry " << iProg << "/" << nEntriesPU << "..." << endl;
      } else {
        if (iProg == nEntriesPU) {
          cout << "        Processing with-pileup entry " << iProg << "/" << nEntriesPU << "..." << endl;
        } else {
          cout << "        Processing with-pileup entry " << iProg << "/" << nEntriesPU << "...\r" << flush;
        }
      }

      // perform calculations
      const Double_t umDcaXY    = pu_dca3dxy * 10000;
      const Double_t umDcaZ     = pu_dca3dz * 10000;
      const Double_t deltaDcaXY = abs(pu_dca3dxysigma / pu_dca3dxy);
      const Double_t deltaDcaZ  = abs(pu_dca3dzsigma / pu_dca3dz);
      const Double_t deltaEta   = abs(pu_deltaeta / pu_eta);
      const Double_t deltaPhi   = abs(pu_deltaphi / pu_phi);
      const Double_t deltaPt    = abs(pu_deltapt / pu_pt);

      // check if values are defined
      const Bool_t thereAreNans = (isnan(pu_dca3dxy) || isnan(pu_dca3dz) || isnan(pu_eta) || isnan(pu_phi) || isnan(pu_pt));
      if (thereAreNans) continue;

      // set reco track variables
      recoTrkVars[TRKVAR::VX]       = pu_vx;
      recoTrkVars[TRKVAR::VY]       = pu_vy;
      recoTrkVars[TRKVAR::VZ]       = pu_vz;
      recoTrkVars[TRKVAR::NMMS]     = (Double_t) pu_nlmms;
      recoTrkVars[TRKVAR::NMAP]     = (Double_t) pu_nlmaps;
      recoTrkVars[TRKVAR::NINT]     = (Double_t) pu_nlintt;
      recoTrkVars[TRKVAR::NTPC]     = (Double_t) pu_ntpc;
      recoTrkVars[TRKVAR::QUAL]     = pu_quality;
      recoTrkVars[TRKVAR::DCAXY]    = umDcaXY;
      recoTrkVars[TRKVAR::DCAZ]     = umDcaZ;
      recoTrkVars[TRKVAR::DELDCAXY] = deltaDcaXY;
      recoTrkVars[TRKVAR::DELDCAZ]  = deltaDcaZ;
      recoTrkVars[TRKVAR::NCLUST]   = 0.;
      recoTrkVars[TRKVAR::AVGCLUST] = 0.;

      // set true track variables
      trueTrkVars[TRKVAR::VX]       = pu_gvx;
      trueTrkVars[TRKVAR::VY]       = pu_gvy;
      trueTrkVars[TRKVAR::VZ]       = pu_gvz;
      trueTrkVars[TRKVAR::NMMS]     = (Double_t) pu_gnlmms;
      trueTrkVars[TRKVAR::NMAP]     = (Double_t) pu_gnlmaps;
      trueTrkVars[TRKVAR::NINT]     = (Double_t) pu_gnlintt;
      trueTrkVars[TRKVAR::NTPC]     = (Double_t) pu_gntpc;
      trueTrkVars[TRKVAR::QUAL]     = pu_quality;
      trueTrkVars[TRKVAR::DCAXY]    = umDcaXY;
      trueTrkVars[TRKVAR::DCAZ]     = umDcaZ;
      trueTrkVars[TRKVAR::DELDCAXY] = deltaDcaXY;
      trueTrkVars[TRKVAR::DELDCAZ]  = deltaDcaZ;
      trueTrkVars[TRKVAR::NCLUST]   = 0.;
      trueTrkVars[TRKVAR::AVGCLUST] = 0.;

      // set reco phys variables
      recoPhysVars[PHYSVAR::PHI]    = pu_phi;
      recoPhysVars[PHYSVAR::ETA]    = pu_eta;
      recoPhysVars[PHYSVAR::PT]     = pu_pt;
      recoPhysVars[PHYSVAR::DELPHI] = deltaPhi;
      recoPhysVars[PHYSVAR::DELETA] = deltaEta;
      recoPhysVars[PHYSVAR::DELPT]  = deltaPt;

      // set true phys variables
      truePhysVars[PHYSVAR::PHI]    = gphi;
      truePhysVars[PHYSVAR::ETA]    = geta;
      truePhysVars[PHYSVAR::PT]     = gpt;
      truePhysVars[PHYSVAR::DELPHI] = deltaPhi;
      truePhysVars[PHYSVAR::DELETA] = deltaEta;
      truePhysVars[PHYSVAR::DELPT]  = deltaPt;

      // check for primary tracks
      const Bool_t isPrimary = (pu_gprimary == 1);

      // fill all histograms
      if (doBeforeCuts) {
        FillTrackHistograms(TYPE::PILEUP, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
        if (isPrimary) {
          FillTrackHistograms(TYPE::PRIMARY, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
        } else {
          FillTrackHistograms(TYPE::NONPRIM, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
        }
      }

      // apply cuts
      const Bool_t isGoodTrk = ApplyCuts(isPrimary, (UInt_t) pu_nlmaps, (UInt_t) pu_ntpc, pu_vz, umDcaXY, umDcaZ, pu_quality);
      if (!isGoodTrk) continue;

      // fill cut histograms
      FillTrackHistograms(TYPE::PILE_CUT, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
      if (isPrimary) {
        FillTrackHistograms(TYPE::PRIM_CUT, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
      } else {
        FillTrackHistograms(TYPE::NONPRIM_CUT, recoTrkVars, trueTrkVars, recoPhysVars, truePhysVars);
      }
    }  // end with-pileup entry loop
    cout << "      Finished with-pileup entry loop." << endl;
  }  // end if (doPileup)

  // normalize histograms if needed
  if (doIntNorm) NormalizeHists();
  return;

}  // end Analyze()



void STrackCutStudy::End() {

  // announce method
  cout << "    Ending:" << endl;

  // set histogram styles
  SetHistStyles();

  // plot labels/directories to save in
  if (makePlots) {
    const TString sLabelAllEO("EmbedOnly_BeforeCuts");
    const TString sLabelCutEO("EmbedOnly_AfterCuts");
    const TString sLabelOddEO("EmbedOnly_WeirdVsNormal");
    const TString sLabelAllPU("WithPileup_BeforeCuts");
    const TString sLabelCutPU("WithPileup_AfterCuts");
    const TString sDirPlotAllEO("EmbedOnlyPlots_BeforeCuts");
    const TString sDirPlotCutEO("EmbedOnlyPlots_AfterCuts");
    const TString sDirPlotOddEO("EmbedOnlyPlots_WeirdVsNormal");
    const TString sDirPlotAllPU("AllWithPileupPlots");
    const TString sDirPlotCutPU("CutWithPileupPlots");

    // track types to plot together
    const Ssiz_t nToDrawAllEO(3);
    const Ssiz_t nToDrawCutEO(3);
    const Ssiz_t nToDrawOddEO(3);
    const Ssiz_t nToDrawAllPU(3);
    const Ssiz_t nToDrawCutPU(3);
    const Int_t  sToDrawAllEO[nToDrawAllEO] = {TYPE::TRACK,    TYPE::TRUTH,     TYPE::WEIRD_ALL};
    const Int_t  sToDrawCutEO[nToDrawCutEO] = {TYPE::TRK_CUT,  TYPE::TRU_CUT,   TYPE::WEIRD_CUT};
    const Int_t  sToDrawOddEO[nToDrawOddEO] = {TYPE::WEIRD_SI_CUT, TYPE::WEIRD_TPC_CUT, TYPE::NORM_CUT};
    const Int_t  sToDrawAllPU[nToDrawAllPU] = {TYPE::PILEUP,   TYPE::PRIMARY,   TYPE::NONPRIM};
    const Int_t  sToDrawCutPU[nToDrawCutPU] = {TYPE::PILE_CUT, TYPE::PRIM_CUT,  TYPE::NONPRIM_CUT};

    // create desired plots
    ConstructPlots(nToDrawCutEO, sToDrawCutEO, sDirPlotCutEO, sLabelCutEO);
    ConstructPlots(nToDrawOddEO, sToDrawOddEO, sDirPlotOddEO, sLabelOddEO);
    if (doBeforeCuts) {
      ConstructPlots(nToDrawAllEO, sToDrawAllEO, sDirPlotAllEO, sLabelAllEO);
    }
    if (doPileup) {
      ConstructPlots(nToDrawCutPU, sToDrawCutPU, sDirPlotCutPU, sLabelCutPU);
      if (doBeforeCuts) ConstructPlots(nToDrawAllPU, sToDrawAllPU, sDirPlotAllPU, sLabelAllPU);
    }
    cout << "      Created plots." << endl;
  }

  // save histograms
  SaveHists();

  // close files
  fOut  -> cd();
  fOut  -> Close();
  fInEO -> cd();
  fInEO -> Close();
  fInPU -> cd();
  fInPU -> Close();
  cout << "      Closed files.\n"
       << "  Finished track cut study!\n"
       << endl;
  return ;

}  // end End()

// end ------------------------------------------------------------------------
