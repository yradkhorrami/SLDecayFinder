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

	m_nSLDecayTotal(0),
	m_nSLDecayBHad(0),
	m_nSLDecayCHad(0),

	m_BHadronIndex{},
	m_CHadronIndex{}
{

    // modify processor description
    _description = "MySLDecayFinder finds semi-leptonic decays inside jets" ;


    // register steering parameters: name, description, class-variable, default value
	registerInputCollection( 	LCIO::MCPARTICLE,
	          			"MCParticleCollection" ,
    	      				"Name of the MCParticle collection"  ,
    	      				m_mcParticleCollection,
    	      				std::string("MCParticle")
    					);

	registerOutputCollection( 	LCIO::MCPARTICLE,
					"SemiLeptonicDecays",
					"Collection of semi-leptonic decays",
					m_SLDecaysCollection,
					std::string("SLDecay")
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

	m_col_SLDecays = new LCCollectionVec(LCIO::MCPARTICLE);


	this->Clear();
	this->ExtractCollections(pLCEvent);
	this->FindSLDecays(pLCEvent);

	m_col_SLDecays->parameters().setValue("nBSLD", (int)m_nSLDecayBHad);
	m_col_SLDecays->parameters().setValue("nCSLD", (int)m_nSLDecayCHad);
	m_col_SLDecays->parameters().setValue("nSLD", (int)m_nSLDecayTotal);
	m_col_SLDecays->parameters().setValues("BHadronIndex", (std::vector<int>)m_BHadronIndex);
	m_col_SLDecays->parameters().setValues("CHadronIndex", (std::vector<int>)m_CHadronIndex);
	pLCEvent->addCollection(m_col_SLDecays, m_SLDecaysCollection);

}



void MySLDecayFinder::check( LCEvent *pLCEvent )
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MySLDecayFinder::end()
{

}

void MySLDecayFinder::Clear()
{

	m_BHadronIndex.clear();
	m_CHadronIndex.clear();
	m_nSLDecayTotal = 0;
	m_nSLDecayBHad = 0;
	m_nSLDecayCHad = 0;
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
	try
	{
		const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_mcParticleCollection);
		for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
		{
			const EVENT::MCParticle *pMCParticle = dynamic_cast<EVENT::MCParticle*>(pLCCollection->getElementAt(i));

			if (NULL == pMCParticle)
				throw EVENT::Exception("Collection type mismatch");

			const int absPdgCode(std::abs(pMCParticle->getPDG()));
			if ((floor(absPdgCode/100)==5) || (floor(absPdgCode/1000)==5))
			{
				for (unsigned int d = 0; d < (pMCParticle->getDaughters()).size(); ++d)
				{
					const int absDauPdgCode(std::abs(((pMCParticle->getDaughters())[d])->getPDG()));
					if ((absDauPdgCode == 12) || (absDauPdgCode == 14) || (absDauPdgCode == 16))
					{
						for (unsigned int o_d = 0; o_d < (pMCParticle->getDaughters()).size(); ++o_d)
						{
							if (std::abs(((pMCParticle->getDaughters())[o_d])->getPDG()) == absDauPdgCode - 1)
							{
								++m_nSLDecayBHad;
								m_BHadronIndex.push_back(i);
							}
						}
					}
				}
			}
			if ((floor(absPdgCode/100)==4) || (floor(absPdgCode/1000)==4))
			{
				for (unsigned int d = 0; d < (pMCParticle->getDaughters()).size(); ++d)
				{
					const int absDauPdgCode(std::abs(((pMCParticle->getDaughters())[d])->getPDG()));
					if ((absDauPdgCode == 12) || (absDauPdgCode == 14) || (absDauPdgCode == 16))
					{
						for (unsigned int o_d = 0; o_d < (pMCParticle->getDaughters()).size(); ++o_d)
						{
							if (std::abs(((pMCParticle->getDaughters())[o_d])->getPDG()) == absDauPdgCode - 1)
							{
								++m_nSLDecayCHad;
								m_CHadronIndex.push_back(i);
							}
						}
					}
				}
			}

		}
		m_nSLDecayTotal = m_nSLDecayBHad + m_nSLDecayCHad;
		streamlog_out(DEBUG) << "Number of Semi-Leptonic decay of B-Hadron: " << m_nSLDecayBHad << std::endl;
		streamlog_out(DEBUG) << "Number of Semi-Leptonic decay of C-Hadron: " << m_nSLDecayCHad << std::endl;
		streamlog_out(DEBUG) << "Total Number of Semi-Leptonic decays: " << m_nSLDecayTotal << std::endl;
	}
	catch (...)
     {
         streamlog_out(WARNING) << "Could not extract Semi-Leptonic decay " << std::endl;
     }
}
