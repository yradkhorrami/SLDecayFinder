#include "MySLDecayFinder.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;


MySLDecayFinder aMySLDecayFinder ;


MySLDecayFinder::MySLDecayFinder() :
	Processor("SLDecayFinder"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),

	m_mcEnergyENu(0.f),
	m_mcEnergyELep(0.f),

	m_nSLDecayTotal(0),
	m_nSLDecayBHad(0),
	m_nSLDecayCHad(0),

	m_Hadron_px(0.f),
	m_Hadron_py(0.f),
	m_Hadron_pz(0.f),
	m_Hadron_E(0.f),
	m_px_vis(0.f),
	m_py_vis(0.f),
	m_pz_vis(0.f),
	m_E_vis(0.f),
	m_SLD_vertex_x(0.f),
	m_SLD_vertex_y(0.f),
	m_SLD_vertex_z(0.f),
	m_EnergyCOM(0.f),
	m_Hadron_mass(0.f),
	m_Mass_vis(0.f),
	m_mcNu_px(0.f),
	m_mcNu_py(0.f),
	m_mcNu_pz(0.f),

	m_pTFile(NULL),
	m_pTTree(NULL)


{

    // modify processor description
    _description = "MySLDecayFinder finds semi-leptonic decays inside jets" ;


    // register steering parameters: name, description, class-variable, default value
	registerInputCollection( LCIO::MCPARTICLE,
	          				"MCParticleCollection" ,
    	      				"Name of the MCParticle collection"  ,
    	      				m_mcParticleCollection,
    	      				std::string("MCParticle")
    						);

	registerProcessorParameter(	"RootFile",
							"Name of the output root file",
							m_rootFile,
							std::string("MySLDecayFinder.root")
							);

}



void MySLDecayFinder::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    // usually a good idea to
    printParameters() ;

	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	this->Clear();

	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
	m_pTTree = new TTree("SLDAnalysisTree", "SLDAnalysisTree");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_pTTree->Branch("EnergyCOM", &m_EnergyCOM, "EnergyCOM/F");
	m_pTTree->Branch("nBHadSLDecay", &m_nSLDecayBHad, "nBHadSLDecay/I");
	m_pTTree->Branch("nCHadSLDecay", &m_nSLDecayCHad, "nCHadSLDecay/I");
	m_pTTree->Branch("nSLDecayTotal", &m_nSLDecayTotal, "nSLDecayTotal/I");
	m_pTTree->Branch("SLD_vertex_x", &m_SLD_vertex_x, "SLD_vertex_x/D");
	m_pTTree->Branch("SLD_vertex_y", &m_SLD_vertex_y, "SLD_vertex_y/D");
	m_pTTree->Branch("SLD_vertex_z", &m_SLD_vertex_z, "SLD_vertex_z/D");
	m_pTTree->Branch("Hadron_E", &m_Hadron_E, "Hadron_E/D");
	m_pTTree->Branch("Hadron_px", &m_Hadron_px, "Hadron_px/D");
	m_pTTree->Branch("Hadron_py", &m_Hadron_py, "Hadron_py/D");
	m_pTTree->Branch("Hadron_pz", &m_Hadron_pz, "Hadron_pz/D");
	m_pTTree->Branch("Hadron_mass", &m_Hadron_mass, "Hadron_mass/D");
	m_pTTree->Branch("E_vis", &m_E_vis, "E_vis/D");
	m_pTTree->Branch("px_vis", &m_px_vis, "px_vis/D");
	m_pTTree->Branch("py_vis", &m_py_vis, "py_vis/D");
	m_pTTree->Branch("pz_vis", &m_pz_vis, "pz_vis/D");
	m_pTTree->Branch("Mass_vis", &m_Mass_vis, "Mass_vis/D");
	m_pTTree->Branch("mcEnergyENu", &m_mcEnergyENu, "mcEnergyENu/D");
	m_pTTree->Branch("mcNu_px", &m_mcNu_px, "mcNu_px/D");
	m_pTTree->Branch("mcNu_py", &m_mcNu_py, "mcNu_py/D");
	m_pTTree->Branch("mcNu_pz", &m_mcNu_pz, "mcNu_pz/D");
	m_pTTree->Branch("mcEnergyELep", &m_mcEnergyELep, "mcEnergyELep/D");
}


void MySLDecayFinder::processRunHeader( LCRunHeader *pLCRunHeader)
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
}



void MySLDecayFinder::processEvent( LCEvent *pLCEvent )
{

	m_nRun = pLCEvent->getRunNumber();
     m_nEvt = pLCEvent->getEventNumber();
     ++m_nEvtSum;

	if ((m_nEvtSum % 100) == 0)
         std::cout << " processed events: " << m_nEvtSum << std::endl;

	this->Clear();
	this->ExtractCollections(pLCEvent);
	this->FindSLDecays(pLCEvent);
	m_pTTree->Fill();

    m_col_SLDecays = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
}



void MySLDecayFinder::check( LCEvent *pLCEvent )
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MySLDecayFinder::end()
{

	m_pTFile->cd();
	m_pTTree->Write();
	m_pTFile->Close();
	delete m_pTFile;

}

void MySLDecayFinder::Clear()
{

	m_EnergyCOM = 0.f;
	m_nSLDecayTotal = 0;
	m_nSLDecayBHad = 0;
	m_nSLDecayCHad = 0;
	m_Hadron_E = 0;
	m_Hadron_px = 0;
	m_Hadron_py = 0;
	m_Hadron_pz = 0;
	m_E_vis = 0;
	m_px_vis = 0;
	m_py_vis = 0;
	m_pz_vis = 0;
	m_SLD_vertex_x = 0;
	m_SLD_vertex_y = 0;
	m_SLD_vertex_z = 0;
	m_mcEnergyELep = 0.f;
	m_mcEnergyENu = 0.f;
	m_Hadron_mass = 0.f;
	m_Mass_vis = 0.f;
	m_mcNu_px = 0.f;
	m_mcNu_py = 0.f;
	m_mcNu_pz = 0.f;

}

void MySLDecayFinder::ExtractCollections(EVENT::LCEvent *pLCEvent)
{
	try
     {
     	const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_mcParticleCollection);

     	for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
		{
			const EVENT::MCParticle *pMCParticle = dynamic_cast<EVENT::MCParticle*>(pLCCollection->getElementAt(i));

			if (NULL == pMCParticle)
			throw EVENT::Exception("Collection type mismatch");

			if (!pMCParticle->getParents().empty())
				continue;

		}
	}
     catch (...)
     {
		streamlog_out(WARNING) << "Could not extract mc particle collection " << m_mcParticleCollection << std::endl;
     }
}

void MySLDecayFinder::FindSLDecays(EVENT::LCEvent *pLCEvent)
{
	MCParticleVector m_mcUnstableParent;
	try
	{
		const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_mcParticleCollection);
		m_EnergyCOM = pLCEvent->getParameters().getFloatVal("Energy");
		for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
		{
			const EVENT::MCParticle *pMCParticle = dynamic_cast<EVENT::MCParticle*>(pLCCollection->getElementAt(i));

			if (NULL == pMCParticle)
				throw EVENT::Exception("Collection type mismatch");

			const int absPdgCode(std::abs(pMCParticle->getPDG()));
			if ((absPdgCode == 12) || (absPdgCode == 14) || (absPdgCode == 16))
			{
				if (!(pMCParticle->getParents().empty()))
				{
					for (unsigned int p = 0; p < (pMCParticle->getParents()).size(); ++p)
					{
						const int absParPdgCode(std::abs((pMCParticle->getParents())[p]->getPDG()));
						if ((floor(absParPdgCode/100)==5) || (floor(absParPdgCode/1000)==5))
						{
							for (unsigned int d = 0; d < (((pMCParticle->getParents())[p])->getDaughters()).size(); ++d)
							{
								const int absDauPdgCode(std::abs(((((pMCParticle->getParents())[p])->getDaughters())[d])->getPDG()));
								if (absDauPdgCode == absPdgCode - 1)
								{
									++m_nSLDecayTotal;
									++m_nSLDecayBHad;
									m_mcEnergyENu += pMCParticle->getEnergy();
									m_mcNu_px += pMCParticle->getMomentum()[0];
									m_mcNu_py += pMCParticle->getMomentum()[1];
									m_mcNu_pz += pMCParticle->getMomentum()[2];
									m_mcEnergyELep += ((((pMCParticle->getParents())[p])->getDaughters())[d])->getEnergy();
									streamlog_out(DEBUG) << "One Semi-Leptonic decay of B-Hadron was found" << std::endl;
									m_mcUnstableParent.push_back((pMCParticle->getParents())[p]);
									streamlog_out(DEBUG) << "found parent of semi-leptonic decay ; nSLD = " << m_mcUnstableParent.size() << std::endl;
									MCParticleVector m_mcDaughtersVector;
									m_SLD_vertex_x = ((pMCParticle->getParents())[p])->getEndpoint()[0];
									m_SLD_vertex_y = ((pMCParticle->getParents())[p])->getEndpoint()[1];
									m_SLD_vertex_z = ((pMCParticle->getParents())[p])->getEndpoint()[2];
									m_Hadron_px = ((pMCParticle->getParents())[p])->getMomentumAtEndpoint()[0];
									m_Hadron_py = ((pMCParticle->getParents())[p])->getMomentumAtEndpoint()[1];
									m_Hadron_pz = ((pMCParticle->getParents())[p])->getMomentumAtEndpoint()[2];
									m_Hadron_E = ((pMCParticle->getParents())[p])->getEnergy();
									m_Hadron_mass = ((pMCParticle->getParents())[p])->getMass();
									for (unsigned int s = 0; s < (((pMCParticle->getParents())[p])->getDaughters()).size(); ++s)
									{
										const EVENT::MCParticle *m_mcDaughters = dynamic_cast<EVENT::MCParticle*>((((pMCParticle->getParents())[p])->getDaughters())[s]);
										if (std::abs(m_mcDaughters->getPDG())!=12 && std::abs(m_mcDaughters->getPDG())!=14 && std::abs(m_mcDaughters->getPDG())!=16)
										{
											m_E_vis += m_mcDaughters->getEnergy();
											m_px_vis += m_mcDaughters->getMomentum()[0];
											m_py_vis += m_mcDaughters->getMomentum()[1];
											m_pz_vis += m_mcDaughters->getMomentum()[2];
											m_Mass_vis += m_mcDaughters->getMass();
										}
									}
								}
							}
						}
						if ((floor(absParPdgCode/100)==4) || (floor(absParPdgCode/1000)==4))
						{
							for (unsigned int d = 0; d < (((pMCParticle->getParents())[p])->getDaughters()).size(); ++d)
							{
								const int absDauPdgCode(std::abs(((((pMCParticle->getParents())[p])->getDaughters())[d])->getPDG()));
								if (absDauPdgCode == absPdgCode - 1)
								{
									++m_nSLDecayTotal;
									++m_nSLDecayCHad;
									m_mcEnergyENu += pMCParticle->getEnergy();
									m_mcNu_px += pMCParticle->getMomentum()[0];
									m_mcNu_py += pMCParticle->getMomentum()[1];
									m_mcNu_pz += pMCParticle->getMomentum()[2];
									m_mcEnergyELep += ((((pMCParticle->getParents())[p])->getDaughters())[d])->getEnergy();
									streamlog_out(DEBUG) << "One Semi-Leptonic decay of Charmed-Hadron was found" << std::endl;
									m_mcUnstableParent.push_back((pMCParticle->getParents())[p]);
									streamlog_out(DEBUG) << "found parent of semi-leptonic decay ; nSLD = " << m_mcUnstableParent.size() << std::endl;
									MCParticleVector m_mcDaughtersVector;
									m_SLD_vertex_x = ((pMCParticle->getParents())[p])->getEndpoint()[0];
									m_SLD_vertex_y = ((pMCParticle->getParents())[p])->getEndpoint()[1];
									m_SLD_vertex_z = ((pMCParticle->getParents())[p])->getEndpoint()[2];
									m_Hadron_px = ((pMCParticle->getParents())[p])->getMomentumAtEndpoint()[0];
									m_Hadron_py = ((pMCParticle->getParents())[p])->getMomentumAtEndpoint()[1];
									m_Hadron_pz = ((pMCParticle->getParents())[p])->getMomentumAtEndpoint()[2];
									m_Hadron_E = ((pMCParticle->getParents())[p])->getEnergy();
									m_Hadron_mass = ((pMCParticle->getParents())[p])->getMass();
									for (unsigned int s = 0; s < (((pMCParticle->getParents())[p])->getDaughters()).size(); ++s)
									{
										const EVENT::MCParticle *m_mcDaughters = dynamic_cast<EVENT::MCParticle*>((((pMCParticle->getParents())[p])->getDaughters())[s]);
										if (std::abs(m_mcDaughters->getPDG())!=12 && std::abs(m_mcDaughters->getPDG())!=14 && std::abs(m_mcDaughters->getPDG())!=16)
										{
											m_E_vis += m_mcDaughters->getEnergy();
											m_px_vis += m_mcDaughters->getMomentum()[0];
											m_py_vis += m_mcDaughters->getMomentum()[1];
											m_pz_vis += m_mcDaughters->getMomentum()[2];
											m_Mass_vis += m_mcDaughters->getMass();
										}
									}
								}
							}
						}
					}
				}
			}
		}
		streamlog_out(DEBUG) << "Number of Semi-Leptonic decay of B-Hadron: " << m_nSLDecayBHad << std::endl;
		streamlog_out(DEBUG) << "Number of Semi-Leptonic decay of C-Hadron: " << m_nSLDecayCHad << std::endl;
		streamlog_out(DEBUG) << "Total Number of Semi-Leptonic decays: " << m_nSLDecayTotal << std::endl;
	}
	catch (...)
     {
         streamlog_out(WARNING) << "Could not extract Semi-Leptonic decay " << std::endl;
     }
}
