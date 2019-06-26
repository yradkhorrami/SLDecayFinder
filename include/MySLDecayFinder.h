#ifndef MySLDecayFinder_h
#define MySLDecayFinder_h 1

#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include "IMPL/LCCollectionVec.h"
#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
//#include "CalibrationHelper.h"

#include <set>
#include <vector>

class TFile;
class TH1F;
class TTree;

using namespace lcio ;
using namespace marlin ;


/**  Example processor for marlin.
 *
 *  If compiled with MARLIN_USE_AIDA
 *  it creates a histogram (cloud) of the MCParticle energies.
 *
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4>
 *  A histogram.
 *
 * @param CollectionName Name of the MCParticle collection
 *
 * @author F. Gaede, DESY
 * @version $Id: MySLDecayFinder.h,v 1.4 2005-10-11 12:57:39 gaede Exp $
 */

class MySLDecayFinder : public Processor {

 public:

  virtual Processor*  newProcessor() { return new MySLDecayFinder ; }


  MySLDecayFinder() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( lcio::LCRunHeader *pLCRunHeader ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( EVENT::LCEvent *pLCEvent ) ;


  virtual void check( EVENT::LCEvent *pLCEvent ) ;

private:

  void Clear();

  void ExtractCollections(EVENT::LCEvent *pLCEvent);

  void FindSLDecays(EVENT::LCEvent *pLCEvent);

//  void CalculateVisible4Momentum();

  /** Called after data processing for clean up.
   */
virtual void end() ;


 protected:

  /** Input collection name.
   */
	std::string 					_colName{} ;

	int						m_nRun{} ;
	int						m_nEvt{} ;
	int						m_nRunSum{};
	int						m_nEvtSum{};

	std::string					m_mcParticleCollection{};
	std::string					m_rootFile{};
	std::string					m_SLDecaysCollection{};
	LCCollectionVec					*m_col_SLDecays{};


	int						m_nSLDecayTotal;
	int						m_nSLDecayBHad;
	int						m_nSLDecayCHad;

	typedef std::vector<int>			IntVector;
	IntVector					m_BHadronIndex;
	IntVector					m_CHadronIndex;
	typedef std::vector<const EVENT::MCParticle*>	MCParticleVector;
	MCParticleVector				m_mcUnstableParent{};
	MCParticleVector				m_mcDaughtersVector{};

} ;

#endif
