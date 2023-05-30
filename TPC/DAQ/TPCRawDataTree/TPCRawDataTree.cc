
#include "TPCRawDataTree.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#include <cassert>
#include <iostream>
#include <memory>

//____________________________________________________________________________..
TPCRawDataTree::TPCRawDataTree(const std::string &name)
  : SubsysReco("TPCRawDataTree")
  , m_fname(name)
{
  // reserve memory for max ADC samples
  m_adcSamples.resize(1024, 0);
}

//____________________________________________________________________________..
int TPCRawDataTree::InitRun(PHCompositeNode *)
{
  m_file = TFile::Open(m_fname.c_str(), "recreate");
  assert(m_file->IsOpen());

  m_PacketTree = new TTree("PacketTree", "Each entry is one packet");

  m_PacketTree->Branch("packet", &m_packet, "packet/I");
  m_PacketTree->Branch("frame", &m_frame, "frame/I");
  m_PacketTree->Branch("nWaveormInFrame", &m_nWaveormInFrame, "nWaveormInFrame/I");
  m_PacketTree->Branch("nTaggerInFrame", &m_nTaggerInFrame, "nTaggerInFrame/I");
  m_PacketTree->Branch("maxFEECount", &m_maxFEECount, "maxFEECount/I");

  m_SampleTree = new TTree("SampleTree", "Each entry is one waveform");

  m_SampleTree->Branch("packet", &m_packet, "packet/I");
  m_SampleTree->Branch("frame", &m_frame, "frame/I");
  m_SampleTree->Branch("nWaveormInFrame", &m_nWaveormInFrame, "nWaveormInFrame/I");
  m_SampleTree->Branch("maxFEECount", &m_maxFEECount, "maxFEECount/I");
  m_SampleTree->Branch("nSamples", &m_nSamples, "nSamples/I");
  m_SampleTree->Branch("adcSamples", &m_adcSamples[0], "adcSamples[nSamples]/s");
  m_SampleTree->Branch("fee", &m_fee, "fee/I");
  m_SampleTree->Branch("sampaAddress", &m_sampaAddress, "sampaAddress/I");
  m_SampleTree->Branch("sampaChannel", &m_sampaChannel, "sampaChannel/I");
  m_SampleTree->Branch("Channel", &m_Channel, "Channel/I");
  m_SampleTree->Branch("BCO", &m_BCO, "BCO/I");
  m_SampleTree->Branch("checksum", &m_checksum, "checksum/I");
  m_SampleTree->Branch("checksumError", &m_checksumError, "checksumError/I");

  m_TaggerTree = new TTree("TaggerTree", "Each entry is one tagger for level 1 trigger or endat tag");

  m_TaggerTree->Branch("packet", &m_packet, "packet/I");
  m_TaggerTree->Branch("frame", &m_frame, "frame/I");
  m_TaggerTree->Branch("tagger_type", &m_tagger_type, "tagger_type/s");
  m_TaggerTree->Branch("is_endat", &m_is_endat, "is_endat/b");
  m_TaggerTree->Branch("is_lvl1", &m_is_lvl1, "is_lvl1/b");
  m_TaggerTree->Branch("bco", &m_bco, "bco/l");
  m_TaggerTree->Branch("lvl1_count", &m_lvl1_count, "lvl1_count/i");
  m_TaggerTree->Branch("endat_count", &m_endat_count, "endat_count/i");
  m_TaggerTree->Branch("last_bco", &m_last_bco, "last_bco/l");
  m_TaggerTree->Branch("modebits", &m_modebits, "modebits/b");

  R1_hist = new TH1F("R1_hist","R1_hist",1024,-0.5,1023.5);
  R2_hist = new TH1F("R2_hist","R2_hist",1024,-0.5,1023.5);
  R3_hist = new TH1F("R3_hist","R3_hist",1024,-0.5,1023.5);

  R1_time = new TH2F("R1_time","R1_time",360,-0.5,359.5,1024,-0.5,1023.5);
  R2_time = new TH2F("R2_time","R2_time",360,-0.5,359.5,1024,-0.5,1023.5);
  R3_time = new TH2F("R3_time","R3_time",360,-0.5,359.5,1024,-0.5,1023.5);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TPCRawDataTree::process_event(PHCompositeNode *topNode)
{
  Event *_event = findNode::getClass<Event>(topNode, "PRDF");
  if (_event == nullptr)
  {
    std::cout << "TPCRawDataTree::Process_Event - Event not found" << std::endl;
    return -1;
  }
  if (_event->getEvtType() >= 8)  /// special events
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  m_frame = _event->getEvtSequence();

  for (int packet : m_packets)
  {
    if (Verbosity())
    {
      std::cout << __PRETTY_FUNCTION__ << " : decoding packet " << packet << std::endl;
    }

    m_packet = packet;

    std::unique_ptr<Packet> p(_event->getPacket(m_packet));
    if (!p)
    {
      if (Verbosity())
      {
        std::cout << __PRETTY_FUNCTION__ << " : missing packet " << packet << std::endl;
      }

      continue;
    }

    m_nWaveormInFrame = p->iValue(0, "NR_WF");
    m_nTaggerInFrame = p->lValue(0, "N_TAGGER");
    m_maxFEECount = p->iValue(0, "MAX_FEECOUNT");

    for (int t = 0; t < m_nTaggerInFrame; t++)
    {
      /*uint16_t*/ m_tagger_type = (uint16_t) (p->lValue(t, "TAGGER_TYPE"));
      /*uint8_t*/ m_is_endat = (uint8_t) (p->lValue(t, "IS_ENDAT"));
      /*uint8_t*/ m_is_lvl1 = (uint8_t) (p->lValue(t, "IS_LEVEL1_TRIGGER"));
      /*uint64_t*/ m_bco = (uint64_t) (p->lValue(t, "BCO"));
      /*uint32_t*/ m_lvl1_count = (uint32_t) (p->lValue(t, "LEVEL1_COUNT"));
      /*uint32_t*/ m_endat_count = (uint32_t) (p->lValue(t, "ENDAT_COUNT"));
      /*uint64_t*/ m_last_bco = (uint64_t) (p->lValue(t, "LAST_BCO"));
      /*uint8_t*/ m_modebits = (uint8_t) (p->lValue(t, "MODEBITS"));

      m_TaggerTree->Fill();
    }

    for (int wf = 0; wf < m_nWaveormInFrame; wf++)
    {
      m_BCO = p->iValue(wf, "BCO");
      m_nSamples = p->iValue(wf, "SAMPLES");
      m_fee = p->iValue(wf, "FEE");

      m_sampaAddress = p->iValue(wf, "SAMPAADDRESS");
      m_sampaChannel = p->iValue(wf, "SAMPACHANNEL");
      m_Channel = p->iValue(wf, "CHANNEL");
      m_checksum = p->iValue(wf, "CHECKSUM");
      m_checksumError = p->iValue(wf, "CHECKSUMERROR");

      TH1F *fillHist;
      TH2F *fillHist2D;

      if(m_fee == 2 ||
         m_fee == 4 ||
         m_fee == 3 ||
         m_fee == 13 ||
         m_fee == 17 ||
         m_fee == 16){ fillHist=R1_hist; fillHist2D=R1_time;}
      else if(m_fee == 11 ||
         m_fee == 12 ||
         m_fee == 19 ||
         m_fee == 18 ||
         m_fee == 01 ||
         m_fee == 00 ||
         m_fee == 16 ||
         m_fee == 15){ fillHist=R2_hist; fillHist2D=R2_time;}
      else{ fillHist=R3_hist; fillHist2D=R3_time;}


      assert(m_nSamples < (int) m_adcSamples.size());  // no need for movements in memory allocation
      for (int s = 0; s < m_nSamples; s++)
      {
        m_adcSamples[s] = p->iValue(wf, s);
        if(m_checksumError==0){
           fillHist->Fill(m_adcSamples[s]);
           fillHist2D->Fill(s,m_adcSamples[s]);
        }
      }

      m_SampleTree->Fill();
    }

    m_PacketTree->Fill();
  }  //   for (int packet : m_packets)

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TPCRawDataTree::End(PHCompositeNode * /*topNode*/)
{
  m_file->Write();

  std::cout << __PRETTY_FUNCTION__ << " : completed saving to " << m_file->GetName() << std::endl;
  m_file->ls();

  m_file->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}
