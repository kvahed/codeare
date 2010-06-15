#include "ReconContext.h"

// Raw strategies 
#include "HannWindow.h"
#include "GenericRegrid.h"
//#include "ConvolutionGridder.h"

// Pixel strategies
#include "MedianFilter.h"
//#include "NUFFT.h"
//#include "CGSENSE.h"

// Raw & Pixel strategies
#include "DummyRecon.h"
#include "InvertOrder.h"
#include "DumpToFile.h"


ReconContext::ReconContext () {

	m_strategies.begin();

	m_strategies.push_back( ((ReconStrategy*) new DummyRecon())   );
	m_strategies.push_back( ((ReconStrategy*) new HannWindow())   );
	m_strategies.push_back( ((ReconStrategy*) new MedianFilter()) );
	m_strategies.push_back( ((ReconStrategy*) new InvertOrder())  );
	m_strategies.push_back( ((ReconStrategy*) new DumpToFile())   );
		
}

