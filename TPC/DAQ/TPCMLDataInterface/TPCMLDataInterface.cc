// $Id: $

/*!
 * \file TPCMLDataInterface.cpp
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "TPCMLDataInterface.h"

#include <trackbase/TpcDefs.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
//#include <g4detectors/PHG4TpcCylinderGeom.h>
//#include <trackbase_historic/SvtxHit.h>
//#include <trackbase_historic/SvtxHitMap.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxEvaluator.h>
#include <g4eval/SvtxHitEval.h>
#include <g4eval/SvtxTruthEval.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

//#include <HepMC/GenEvent.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <H5Cpp.h>

#include <boost/format.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;
using namespace CLHEP;
using namespace H5;

TPCMLDataInterface::TPCMLDataInterface(
    unsigned int minLayer,
    unsigned int m_maxLayer,
    const std::string& outputfilename)
  : SubsysReco("TPCMLDataInterface")
  , m_svtxevalstack(nullptr)
  , m_strict(false)
  , m_saveDataStreamFile(true)
  , m_outputFileNameBase(outputfilename)
  , m_h5File(nullptr)
  , m_minLayer(minLayer)
  , m_maxLayer(m_maxLayer)
  , m_evtCounter(-1)
  , m_vertexZAcceptanceCut(10)
  , m_etaAcceptanceCut(1.1)
  , m_momentumCut(.1)
  , m_hDataSize(nullptr)
  , m_hWavelet(nullptr)
  , m_hNChEta(nullptr)
  , m_hLayerWaveletSize(nullptr)
  , m_hLayerHit(nullptr)
  , m_hLayerZBinHit(nullptr)
  , m_hLayerZBinADC(nullptr)
  , m_hLayerDataSize(nullptr)
  , m_hLayerSumHit(nullptr)
  , m_hLayerSumDataSize(nullptr)
{
}

TPCMLDataInterface::~TPCMLDataInterface()
{
  if (m_h5File) delete m_h5File;
}

int TPCMLDataInterface::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCMLDataInterface::End(PHCompositeNode* topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
    cout << "TPCMLDataInterface::End - write to " << m_outputFileNameBase + ".root" << endl;
  PHTFileServer::get().cd(m_outputFileNameBase + ".root");
  /*
  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);
  for (unsigned int i = 0; i < hm->nHistos(); i++)
    hm->getHisto(i)->Write();

  // help index files with TChain
  TTree* T_Index = new TTree("T_Index", "T_Index");
  assert(T_Index);
  T_Index->Write();
  */
  _rawHits->Write();
  h_TXY->Write();
  /*
  if (m_h5File)
  {
    delete m_h5File;
    m_h5File = nullptr;
  }
  */
  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCMLDataInterface::InitRun(PHCompositeNode* topNode)
{
  PHG4TpcCylinderGeomContainer* seggeo = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!seggeo)
  {
    cout << "could not locate geo node "
         << "CYLINDERCELLGEOM_SVTX" << endl;
    exit(1);
  }

  int nZBins = 0;
  for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    PHG4TpcCylinderGeom* layerGeom =
        seggeo->GetLayerCellGeom(layer);
    assert(layerGeom);

    if (nZBins <= 0)
    {
      nZBins = layerGeom->get_zbins();
      assert(nZBins > 0);
    }
    else
    {
      if ((int) nZBins != layerGeom->get_zbins())
      {
        cout << "TPCMLDataInterface::InitRun - Fatal Error - nZBin at layer " << layer << " is " << layerGeom->get_zbins()
             << ", which is different from previous layers of nZBin = " << nZBins << endl;
        exit(1);
      }
    }
  }  //   for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  if (Verbosity() >= VERBOSITY_SOME)
    cout << "TPCMLDataInterface::get_HistoManager - Making PHTFileServer " << m_outputFileNameBase + ".root"
         << endl;
  PHTFileServer::get().open(m_outputFileNameBase + ".root", "RECREATE");

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);

  TH1D* h = new TH1D("hNormalization",  //
                     "Normalization;Items;Summed quantity", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event count");
  h->GetXaxis()->SetBinLabel(i++, "Collision count");
  h->GetXaxis()->SetBinLabel(i++, "TPC G4Hit");
  h->GetXaxis()->SetBinLabel(i++, "TPC G4Hit Edep");
  h->GetXaxis()->SetBinLabel(i++, "TPC Hit");
  h->GetXaxis()->SetBinLabel(i++, "TPC Wavelet");
  h->GetXaxis()->SetBinLabel(i++, "TPC DataSize");

  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  //  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  //  {
  //    const PHG4TpcCylinderGeom* layer_geom = seggeo->GetLayerCellGeom(layer);

  //    const string histNameCellHit(boost::str(boost::format{"hCellHit_Layer%1%"} % layer));
  //    const string histNameCellCharge(boost::str(boost::format{"hCellCharge_Layer%1%"} % layer));

  //  }

  hm->registerHisto(m_hDataSize =
                        new TH1D("hDataSize",  //
                                 "TPC Data Size per Event;Data size [Byte];Count",
                                 10000, 0, 20e6));

  hm->registerHisto(m_hWavelet =
                        new TH1D("hWavelet",  //
                                 "TPC Recorded Wavelet per Event;Wavelet count;Count",
                                 10000, 0, 4e6));

  hm->registerHisto(m_hNChEta =
                        new TH1D("hNChEta",  //
                                 "Charged particle #eta distribution;#eta;Count",
                                 1000, -5, 5));

  hm->registerHisto(m_hLayerWaveletSize =
                        new TH2D("hLayerWaveletSize",  //
                                 "Number of Recorded ADC sample per Wavelet;Layer ID;ADC Sample Count per Wavelet",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 nZBins, -.5, nZBins - .5));

  hm->registerHisto(m_hLayerHit =
                        new TH2D("hLayerHit",  //
                                 "Number of Recorded ADC sample per channel;Layer ID;ADC Sample Count",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 nZBins, -.5, nZBins - .5));

  hm->registerHisto(m_hLayerDataSize =
                        new TH2D("hLayerDataSize",  //
                                 "Data size per channel;Layer ID;Data size [Byte]",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 2 * nZBins, 0, 2 * nZBins));

  hm->registerHisto(m_hLayerSumHit =
                        new TH2D("hLayerSumHit",  //
                                 "Number of Recorded ADC sample per layer;Layer ID;ADC Sample Count",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 10000, -.5, 99999.5));

  hm->registerHisto(m_hLayerSumDataSize =
                        new TH2D("hLayerSumDataSize",  //
                                 "Data size per trigger per layer;Layer ID;Data size [Byte]",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 1000, 0, .5e6));

  hm->registerHisto(m_hLayerZBinHit =
                        new TH2D("hLayerZBinHit",  //
                                 "Number of Recorded ADC sample per Time Bin;z bin ID;Layer ID",
                                 nZBins, -.5, nZBins - .5,
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5));

  hm->registerHisto(m_hLayerZBinADC =
                        new TH2D("hLayerZBinADC",  //
                                 "Sum ADC per Time Bin;z bin ID;Layer ID",
                                 nZBins, -.5, nZBins - .5,
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5));

  // init HDF5 output

  if (Verbosity() >= VERBOSITY_SOME)
    cout << "TPCMLDataInterface::get_HistoManager - Making H5File " << m_outputFileNameBase + ".h5"
         << endl;
  //m_h5File = new H5File(m_outputFileNameBase + ".h5", H5F_ACC_TRUNC);
  // side sector time R Pad
  //float adc_R1[2][12][260][16][1152] ;
  //float adc_R2[2][12][260][16][1536] ;
  //float adc_R3[2][12][260][16][2304] ;
  memset( adc_R1, 0, sizeof(adc_R1) ); // initialize all the elments to zero
  memset( adc_R2, 0, sizeof(adc_R2) ); // initialize all the elments to zero
  memset( adc_R3, 0, sizeof(adc_R3) ); // initialize all the elments to zero

  //Int_t nsides = 2;
  //Int_t nsectors = 24;
  //Int_t ntimebins = 260;
  //Int_t nlayers = 16;
  //Int_t npadsr1 = 96;
  //Int_t npadsr2 = 128;
  //Int_t npadsr3 = 192;

  _hit_r = 0.0;
  _hit_adc  = 0;
  // TTree
  _rawHits = new TTree("RawHitsTree", "tpc hit tree for ionization");
  //_rawHits->Branch("nsides",   &nsides   ,"nsides/I");

  _rawHits->Branch("adc_R1", adc_R1,"adc_R1[24][260][16][96]/I");
  _rawHits->Branch("adc_R2", adc_R2,"adc_R2[24][260][16][128]/I");
  _rawHits->Branch("adc_R3", adc_R3,"adc_R3[24][260][16][192]/I");
  h_TXY = new TH3F("h_TXY","h_TXY",260,-0.5,259.5,120,-60,60,120,-60,60);

  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCMLDataInterface::process_event(PHCompositeNode* topNode)
{
  m_evtCounter += 1;

  //assert(m_h5File);

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);
  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto("hNormalization"));
  assert(h_norm);
  h_norm->Fill("Event count", 1);

  if (!m_svtxevalstack)
  {
    m_svtxevalstack = new SvtxEvalStack(topNode);
    m_svtxevalstack->set_strict(m_strict);
    m_svtxevalstack->set_verbosity(Verbosity());
    m_svtxevalstack->set_use_initial_vertex(m_use_initial_vertex);
    m_svtxevalstack->set_use_genfit_vertex(m_use_genfit_vertex);
    m_svtxevalstack->next_event(topNode);
  }
  else
  {
    m_svtxevalstack->next_event(topNode);
  }

  //SvtxHitEval* hiteval = m_svtxevalstack->get_hit_eval();
  //SvtxTruthEval* trutheval = m_svtxevalstack->get_truth_eval();

  //float gpx = 0;
  //float gpy = 0;
  //float gpz = 0;

  PHG4HitContainer* g4hit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  cout << "TPCMLDataInterface::process_event - g4 hit node G4HIT_TPC" << endl;
  if (!g4hit)
  {
    cout << "TPCMLDataInterface::process_event - Could not locate g4 hit node G4HIT_TPC" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node containing the digitized hits
  TrkrHitSetContainer* hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hits)
  {
    cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  //TrkrHitSetContainerv1* hits_v1 = findNode::getClass<TrkrHitSetContainerv1>(topNode, "TRKR_HITSET");
  //if (!hits_v1)
  //{
  //  cout << PHWHERE << "TrkrHitSetContainerv1 ERROR: Can't find node TRKR_HITSET" << endl;
  //  return Fun4AllReturnCodes::ABORTRUN;
  //}

  //  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_SVTX");
  //  if (!cells)
  //  {
  //    cout << "TPCMLDataInterface::process_event - could not locate cell node "
  //         << "G4CELL_SVTX" << endl;
  //    exit(1);
  //  }

  PHG4TpcCylinderGeomContainer* seggeo = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!seggeo)
  {
    cout << "TPCMLDataInterface::process_event - could not locate geo node "
         << "CYLINDERCELLGEOM_SVTX" << endl;
    exit(1);
  }

  PHHepMCGenEventMap* geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    static bool once = true;
    if (once)
    {
      once = false;

      cout << "TPCMLDataInterface::process_event - - missing node PHHepMCGenEventMap. Skipping HepMC stat." << std::endl;
    }
  }
  else
  {
    h_norm->Fill("Collision count", geneventmap->size());
  }

  PHG4TruthInfoContainer* truthInfoList = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                                     "G4TruthInfo");
  if (!truthInfoList)
  {
    cout << "TPCMLDataInterface::process_event - Fatal Error - "
         << "unable to find DST node "
         << "G4TruthInfo" << endl;
    assert(truthInfoList);
  }

  PHG4TruthInfoContainer::ConstRange primary_range =
      truthInfoList->GetPrimaryParticleRange();

  for (PHG4TruthInfoContainer::ConstIterator particle_iter =
           primary_range.first;
       particle_iter != primary_range.second;
       ++particle_iter)
  {
    const PHG4Particle* p = particle_iter->second;
    assert(p);

    TParticlePDG* pdg_p = TDatabasePDG::Instance()->GetParticle(
        p->get_pid());
    assert(pdg_p);

    if (fabs(pdg_p->Charge()) > 0)
    {
      TVector3 pvec(p->get_px(), p->get_py(), p->get_pz());

      if (pvec.Perp2() > 0)
      {
        assert(m_hNChEta);
        m_hNChEta->Fill(pvec.PseudoRapidity());
      }
    }

  }  //          if (_load_all_particle) else

  for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(layer);

    for (PHG4HitContainer::ConstIterator hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      const double edep = hiter->second->get_edep();
      h_norm->Fill("TPC G4Hit Edep", edep);
      h_norm->Fill("TPC G4Hit", 1);
    }  //     for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)

  }  //   for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  //H5 init
  string h5GroupName = boost::str(boost::format("TimeFrame_%1%") % m_evtCounter);
  //unique_ptr<Group> h5Group(new Group(m_h5File->createGroup("/" + h5GroupName)));
  //h5Group->setComment(boost::str(boost::format("Collection of ADC data matrixes in Time Frame #%1%") % m_evtCounter));
  //if (Verbosity())
  //  cout << "TPCMLDataInterface::process_event - save to H5 group " << h5GroupName << endl;
  map<int, shared_ptr<DataSet>> layerH5DataSetMap;
  map<int, shared_ptr<DataSet>> layerH5SignalBackgroundMap;
  map<int, shared_ptr<DataSpace>> layerH5DataSpaceMap;
  map<int, shared_ptr<DataSpace>> layerH5SignalBackgroundDataSpaceMap;

  //wavelet size stat.
  vector<uint32_t> layerWaveletDataSize(m_maxLayer + 1, 0);
  vector<hsize_t> layerSize(1);
  layerSize[0] = static_cast<hsize_t>(m_maxLayer - m_minLayer + 1);
  shared_ptr<DataSpace> H5DataSpaceLayerWaveletDataSize(new DataSpace(1, layerSize.data()));
  //shared_ptr<DataSet> H5DataSetLayerWaveletDataSize(new DataSet(h5Group->createDataSet(
  //    "sPHENIXRawDataSizeBytePerLayer",
  //    PredType::NATIVE_UINT32,
  //    *(H5DataSpaceLayerWaveletDataSize))));

  // prepreare stat. storage
  int nZBins = 0;
  vector<array<vector<int>, 2>> layerChanHit(m_maxLayer + 1);
  vector<array<vector<int>, 2>> layerChanDataSize(m_maxLayer + 1);
  vector<vector<uint16_t>> layerDataBuffer(m_maxLayer + 1);
  vector<vector<uint8_t>> layerSignalBackgroundBuffer(m_maxLayer + 1);
  vector<vector<hsize_t>> layerDataBufferSize(m_maxLayer + 1, vector<hsize_t>({0, 0}));
  vector<vector<hsize_t>> layerSignalBackgroundDataBufferSize(m_maxLayer + 1, vector<hsize_t>({0, 0}));

  //int nWavelet = 0;
  //int sumDataSize = 0;
  for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    PHG4TpcCylinderGeom* layerGeom =
        seggeo->GetLayerCellGeom(layer);
    assert(layerGeom);

    // start with an empty vector of cells for each phibin
    const int nphibins = layerGeom->get_phibins();
    assert(nphibins > 0);

    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "TPCMLDataInterface::process_event - init layer " << layer << " with "
           << "nphibins = " << nphibins
           << ", layerGeom->get_zbins() = " << layerGeom->get_zbins() << endl;
    }

    if (nZBins <= 0)
    {
      nZBins = layerGeom->get_zbins();
      assert(nZBins > 0);
    }
    else
    {
      if ((int) nZBins != layerGeom->get_zbins())
      {
        cout << "TPCMLDataInterface::process_event - Fatal Error - nZBin at layer " << layer << " is " << layerGeom->get_zbins()
             << ", which is different from previous layers of nZBin = " << nZBins << endl;
        exit(1);
      }
    }

    for (unsigned int side = 0; side < 2; ++side)
    {
      layerChanHit[layer][side].resize(nphibins, 0);

      layerChanDataSize[layer][side].resize(nphibins, 0);
    }  //     for (unsigned int side = 0; side < 2; ++side)

    // buffer init
    layerDataBufferSize[layer][0] = nphibins;
    layerDataBufferSize[layer][1] = layerGeom->get_zbins();
    layerSignalBackgroundDataBufferSize[layer][0] = nphibins;
    layerSignalBackgroundDataBufferSize[layer][1] = layerGeom->get_zbins();
    layerDataBuffer[layer].resize(layerDataBufferSize[layer][0] * layerDataBufferSize[layer][1], 0);
    layerSignalBackgroundBuffer[layer].resize(layerSignalBackgroundDataBufferSize[layer][0] * layerSignalBackgroundDataBufferSize[layer][1], 0);

    static const vector<hsize_t> cdims({32, 32});
    DSetCreatPropList ds_creatplist;
    ds_creatplist.setChunk(2, cdims.data());  // then modify it for compression
    ds_creatplist.setDeflate(6);

    layerH5DataSpaceMap[layer] = shared_ptr<DataSpace>(new DataSpace(2, layerDataBufferSize[layer].data()));
    layerH5SignalBackgroundDataSpaceMap[layer] = shared_ptr<DataSpace>(new DataSpace(2, layerSignalBackgroundDataBufferSize[layer].data()));

    //layerH5DataSetMap[layer] = shared_ptr<DataSet>(new DataSet(h5Group->createDataSet(
    //    boost::str(boost::format("Data_Layer%1%") % (layer - m_minLayer)),
    //    PredType::NATIVE_UINT16,
    //    *(layerH5DataSpaceMap[layer]),
    //    ds_creatplist)));

    //layerH5SignalBackgroundMap[layer] = shared_ptr<DataSet>(new DataSet(h5Group->createDataSet(
    //    boost::str(boost::format("Data_Layer_SignalBackground%1%") % (layer - m_minLayer)),
    //    PredType::NATIVE_UINT8,
    //    *(layerH5SignalBackgroundDataSpaceMap[layer]),
    //    ds_creatplist)));

  }  //   for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  assert(nZBins > 0);

  // count hits and make wavelets
  //int last_layer = -1;
  //int last_side = -1;
  //int last_phibin = -1;
  //int last_zbin = -1;
  vector<unsigned int> last_wavelet;
  //int last_wavelet_hittime = -1;

  //  for (SvtxHitMap::Iter iter = hits->begin(); iter != hits->end(); ++iter)
  //  {
  //    SvtxHit* hit = iter->second;
  // loop over the TPC HitSet objects
  //TrkrHitSetContainer::ConstRange hitsetrange = hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  //eshulga new:
  int ntpc_phibins[3]= {1152/12, 1536/12, 2304/12}; 
  for(int n_sec=0; n_sec<3; n_sec++){
    unsigned int pads_per_sector = ntpc_phibins[n_sec];
    unsigned int sector = 200 / pads_per_sector;
    unsigned int layernum = 1;
    unsigned int side = 1;

    TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layernum, sector, side);
    TrkrHitSetContainer::Iterator hitsetit = hits->findOrAddHitSet(hitsetkey);

    const int layer_tmp = TrkrDefs::getLayer(hitsetit->first);
    cout << "layer_tmp = " << layer_tmp << endl;
  } 
  for(int side=0;side<2;side++){
    for(int sector=0;sector<12;sector++){
      //for(int layernum=7;layernum<=7+16*3;layernum++){
      for(int layernum=m_minLayer;layernum<=m_maxLayer;layernum++){
        PHG4TpcCylinderGeom* layerGeom =
              seggeo->GetLayerCellGeom(layernum);
        TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layernum, sector, side);
        TrkrHitSetContainer::Iterator hitsetit = hits->findOrAddHitSet(hitsetkey);
        TrkrHitSet* hitset = hitsetit->second;
        TrkrHitSet::ConstRange hitrangei = hitset->getHits();
        for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
             hitr != hitrangei.second;
             ++hitr)
        {

          //int phibin = TpcDefs::getPad(hitr->first);
          int zbin = TpcDefs::getTBin(hitr->first);
          if(zbin>260) continue;
          const double z_abs = fabs(layerGeom->get_zcenter(zbin));
          //cout << "ZBin: " << zbin << " t=" << z_abs << endl;
          const double r = layerGeom->get_radius();
          TVector3 acceptanceVec(r, 0, z_abs - m_vertexZAcceptanceCut);
          //const double eta = acceptanceVec.PseudoRapidity();
          // frame of 249 z bins [r] 
          _hit_r = layerGeom->get_radius();
          _hit_adc = hitr->second->getAdc();
          int phibin = TpcDefs::getPad(hitr->first);
          //if(_hit_adc>0) cout<<"t="<<zbin<<" ADC="<<_hit_adc<<endl;
          //cout<<"side: "<<side<<" sector: "<< sector <<" zbin:"<<zbin<< " layernum: " << layernum <<" layer: " << (layernum-7)/16<<" phibin:"<< phibin <<endl;
          if(zbin<260){
            int Rsec = -1;
            if(layernum<7+16){
              //adc_R1[side][sector][zbin][(layernum-7)][phibin_r1]=_hit_adc;
              Rsec = 0;
              }                        
              
            if(layernum>7+16-1 && layernum<7+2*16)   {
              //adc_R2[side][sector][zbin][(layernum-7)-16][phibin_r2]=_hit_adc;
              Rsec = 1;
              }
            if(layernum>7+2*16-1 && layernum<7+3*16) {
              //adc_R3[side][sector][zbin][(layernum-7)-2*16][phibin_r3]=_hit_adc;
              Rsec = 2;
              }
            int phibin_r = phibin-sector*ntpc_phibins[Rsec];
            //int phibin_r2 = phibin-sector*ntpc_phibins[1];
            //int phibin_r3 = phibin-sector*ntpc_phibins[2];

            if(Rsec==0)adc_R1[11*side+sector][zbin][(layernum-7)-Rsec*16][phibin_r]=_hit_adc;
            if(Rsec==1)adc_R2[11*side+sector][zbin][(layernum-7)-Rsec*16][phibin_r]=_hit_adc;
            if(Rsec==2)adc_R3[11*side+sector][zbin][(layernum-7)-Rsec*16][phibin_r]=_hit_adc;
            h_TXY->Fill(zbin,layernum*cos(phibin*M_PI/(ntpc_phibins[Rsec]*6)),layernum*sin(phibin*M_PI/(ntpc_phibins[Rsec]*6)));
            //cout<<"sector="<<sector<<" phibin"<<  phibin
            //<<"phibin_r1="<<phibin_r1
            //<<"phibin_r2="<<phibin_r2
            //<<"phibin_r3="<<phibin_r3
            //<<endl;
            //cout<<"adc_R1="<<adc_R1[side][sector][zbin][(layernum-7)/16][phibin]<<"_hit_adc"<<_hit_adc<<endl;
            //cout<<"adc_R2="<<adc_R2[side][sector][zbin][(layernum-7)/16][phibin]<<"_hit_adc"<<_hit_adc<<endl;
            //cout<<"adc_R3="<<adc_R3[side][sector][zbin][(layernum-7)/16][phibin]<<"_hit_adc"<<_hit_adc<<endl;
          }
          //if(number % 249 == 0)
                    

        }
      }      
    }
  }
  _rawHits->Fill();
  memset( adc_R1, 0, sizeof(adc_R1) ); // initialize all the elments to zero
  memset( adc_R2, 0, sizeof(adc_R2) ); // initialize all the elments to zero
  memset( adc_R3, 0, sizeof(adc_R3) ); // initialize all the elments to zero
  /*
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    TrkrHitSet* hitset = hitsetitr->second;

    if (Verbosity() > 2) hitset->identify();

    //    const int layer = hit->get_layer();
    // we have a single hitset, get the info that identifies the module
    const int layer = TrkrDefs::getLayer(hitsetitr->first);

    if (Verbosity() >= 2)
    {
      cout << "TPCMLDataInterface::process_event - TrkrHitSet @ layer " << layer << " with " << hitset->size() << " hits" << endl;
    }

    if (layer < m_minLayer or layer > m_maxLayer) continue;

    //    PHG4Cell* cell = cells->findCell(hit->get_cellid());                           //not needed once geofixed

    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
    {
      //    const int phibin = PHG4CellDefs::SizeBinning::get_phibin(cell->get_cellid());  //cell->get_binphi();
      //    const int zbin = PHG4CellDefs::SizeBinning::get_zbin(cell->get_cellid());      //cell->get_binz();
      int phibin = TpcDefs::getPad(hitr->first);
      int zbin = TpcDefs::getTBin(hitr->first);
      const int side = (zbin < nZBins / 2) ? 0 : 1;

      if (Verbosity() >= 2)
      {
        cout << "TPCMLDataInterface::process_event - hit @ layer " << layer << " phibin = " << phibin << " zbin = " << zbin << " side = " << side << endl;
      }

      assert(phibin >= 0);
      assert(zbin >= 0);
      assert(side >= 0);

      // new wavelet?
      if (last_layer != layer or last_phibin != phibin or last_side != side or abs(last_zbin - zbin) != 1)
      {
        // save last wavelet
        if (last_wavelet.size() > 0)
        {
          const int datasize = writeWavelet(last_layer, last_side, last_phibin, last_wavelet_hittime, last_wavelet);
          assert(datasize > 0);

          nWavelet += 1;
          sumDataSize += datasize;
          layerChanDataSize[last_layer][last_side][last_phibin] += datasize;

          last_wavelet.clear();
          last_zbin = -1;
        }

        // z-R cut on digitized wavelet
        PHG4TpcCylinderGeom* layerGeom =
            seggeo->GetLayerCellGeom(layer);
        assert(layerGeom);

        const double z_abs = fabs(layerGeom->get_zcenter(zbin));
        const double r = layerGeom->get_radius();
        TVector3 acceptanceVec(r, 0, z_abs - m_vertexZAcceptanceCut);
        const double eta = acceptanceVec.PseudoRapidity();

        //_hit_r = layerGeom->get_radius();
        //_hit_adc = hitr->second->getAdc();
        //_rawHits->Fill();
        // Event
        //  - side
        //    - sector (0-11)
        //      -> TTree Fill
        //Sector by sector R phi z (249 z bins)
        if (eta > m_etaAcceptanceCut) continue;

        // make new wavelet
        last_layer = layer;
        last_side = side;
        last_phibin = phibin;

        // time check
        last_wavelet_hittime = (side == 0) ? (zbin) : (nZBins - 1 - zbin);
        assert(last_wavelet_hittime >= 0);
        assert(last_wavelet_hittime <= nZBins / 2);

      }  //     if (last_layer != layer or last_phibin != phibin)

      if (Verbosity() >= VERBOSITY_A_LOT)
      {
        cout << "TPCMLDataInterface::process_event -  layer " << layer << " hit with "

             << "phibin = " << phibin
             << ",zbin = " << zbin
             << ",side = " << side
             << ",last_wavelet.size() = " << last_wavelet.size()
             << ",last_zbin = " << last_zbin
             << endl;
      }

      // more checks on signal continuity
      if (last_wavelet.size() > 0)
      {
        //        if (side == 0)
        //        {
        assert(zbin - last_zbin == 1);
        //        }
        //        else
        //        {
        //          assert(last_zbin - zbin == 1);
        //        }
      }

      // record adc
      unsigned int adc = hitr->second->getAdc();
      last_wavelet.push_back(adc);
      last_zbin = zbin;


      // statistics
      layerChanHit[layer][side][phibin] += 1;
      assert(m_hLayerZBinHit);
      m_hLayerZBinHit->Fill(zbin, layer, 1);
      assert(m_hLayerZBinADC);
      m_hLayerZBinADC->Fill(zbin, layer, adc);
      assert(adc < (1 << 10));
      assert((hsize_t) phibin < layerDataBufferSize[layer][0]);
      assert((hsize_t) zbin < layerDataBufferSize[layer][1]);
      const size_t hitindex(layerDataBufferSize[layer][1] * phibin + zbin);
      assert((hsize_t) phibin < layerSignalBackgroundDataBufferSize[layer][0]);
      assert((hsize_t) zbin < layerSignalBackgroundDataBufferSize[layer][1]);
      const size_t hitindexSB(layerSignalBackgroundDataBufferSize[layer][1] * phibin + zbin);
      assert(hitindex < layerDataBuffer[layer].size());
      assert(hitindexSB < layerSignalBackgroundBuffer[layer].size());

      if (layerDataBuffer[layer][hitindex] != 0)
      {
        cout << "TPCMLDataInterface::process_event - WARNING - hit @ layer "
             << layer << " phibin = " << phibin << " zbin = " << zbin << " side = " << side
             << " overwriting previous hit with ADC = " << layerDataBuffer[layer][hitindex]
             << endl;
      }
      layerDataBuffer[layer][hitindex] = adc;

      if (layerSignalBackgroundBuffer[layer][hitindexSB] != 0)
      {
        cout << "TPCMLDataInterface::process_event - WARNING - signal/background @ layer "
             << layer << " phibin = " << phibin << " zbin = " << zbin << " side = " << side
             << " overwriting previous hit with = " << layerSignalBackgroundBuffer[layer][hitindexSB]
             << endl;
      }

      // momentum cut

      bool signal_or_bgd = false;

      TrkrDefs::hitkey hit_key = hitr->first;
      //        PHG4Hit* g4hit = hiteval->max_truth_hit_by_energy(hit_key);
      std::set<PHG4Hit*> g4hit_set = hiteval->all_truth_hits(hit_key);

      for (PHG4Hit* g4hit : g4hit_set)
      {
        PHG4Particle* g4particle = trutheval->get_particle(g4hit);

        if (g4particle != nullptr)
        {
          gpx = g4particle->get_px();
          gpy = g4particle->get_py();
          gpz = g4particle->get_pz();
        }

        if (sqrt(gpx * gpx + gpy * gpy + gpz * gpz) > m_momentumCut) signal_or_bgd = true;

      }  //        for (auto g4hit: g4hit_set)

      layerSignalBackgroundBuffer[layer][hitindexSB] = signal_or_bgd;

    }  //     for (TrkrHitSet::ConstIterator hitr = hitrangei.first;

  }  //   for(SvtxHitMap::Iter iter = hits->begin(); iter != hits->end(); ++iter) {
  
  // save last wavelet
  if (last_wavelet.size() > 0)
  {
    const int datasize = writeWavelet(last_layer, last_side, last_phibin, last_wavelet_hittime, last_wavelet);
    assert(datasize > 0);

    nWavelet += 1;
    sumDataSize += datasize;
    layerChanDataSize[last_layer][last_side][last_phibin] += datasize;
  }

  // statistics
  for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    for (unsigned int side = 0; side < 2; ++side)
    {
      int sumHit = 0;
      for (const int& hit : layerChanHit[layer][side])
      {
        sumHit += hit;

        assert(m_hLayerHit);
        m_hLayerHit->Fill(layer, hit);
        h_norm->Fill("TPC Hit", hit);

        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << "TPCMLDataInterface::process_event - layerChanHit: "
               << "hit = " << hit
               << "at layer = " << layer
               << ", side = " << side
               << endl;
        }
      }

      if (Verbosity() >= VERBOSITY_MORE)
      {
        cout << "TPCMLDataInterface::process_event - hLayerSumCellHit->Fill(" << layer << ", " << sumHit << ")" << endl;
      }
      assert(m_hLayerSumHit);
      m_hLayerSumHit->Fill(layer, sumHit);

      double sumData = 0;
      for (const int& data : layerChanDataSize[layer][side])
      {
        sumData += data;

        assert(m_hLayerDataSize);
        m_hLayerDataSize->Fill(layer, data);
      }

      assert(m_hLayerSumDataSize);
      m_hLayerSumDataSize->Fill(layer, sumData);
      layerWaveletDataSize[(layer - m_minLayer)] += sumData;
    }  //    for (unsigned int side = 0; side < 2; ++side)

    // store in H5

    assert(layerH5DataSetMap[layer]);
    assert(layerH5SignalBackgroundMap[layer]);
    layerH5DataSetMap[layer]->write(layerDataBuffer[layer].data(), PredType::NATIVE_UINT16);
    layerH5SignalBackgroundMap[layer]->write(layerSignalBackgroundBuffer[layer].data(), PredType::NATIVE_UINT8);
    H5DataSetLayerWaveletDataSize->write(layerWaveletDataSize.data(), PredType::NATIVE_UINT32);

  }  //  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  assert(m_hWavelet);
  m_hWavelet->Fill(nWavelet);
  h_norm->Fill("TPC Wavelet", nWavelet);
  assert(m_hDataSize);
  m_hDataSize->Fill(sumDataSize);
  h_norm->Fill("TPC DataSize", sumDataSize);
  */
  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCMLDataInterface::writeWavelet(int layer, int side, int phibin, int hittime, const vector<unsigned int>& wavelet)
{
  static const int headersize = 2;  // 2-byte header per wavelet

  //data in byte aligned and padding
  const int datasizebit = wavelet.size() * 10;
  int datasizebyte = datasizebit / 8;
  if (datasizebyte * 8 < datasizebit) datasizebyte += 1;

  assert(m_hLayerWaveletSize);
  m_hLayerWaveletSize->Fill(layer, wavelet.size());

  return headersize + datasizebyte;
}

Fun4AllHistoManager*
TPCMLDataInterface::getHistoManager()
{
  static string histname("TPCMLDataInterface_HISTOS");

  Fun4AllServer* se = Fun4AllServer::instance();
  Fun4AllHistoManager* hm = se->getHistoManager(histname);

  if (not hm)
  {
    cout
        << "TPCMLDataInterface::get_HistoManager - Making Fun4AllHistoManager " << histname
        << endl;
    hm = new Fun4AllHistoManager(histname);
    se->registerHistoManager(hm);
  }

  assert(hm);

  return hm;
}
