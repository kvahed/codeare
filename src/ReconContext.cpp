/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

#include "ReconContext.h"

// Raw strategies 
#include "HannWindow.h"
#include "GenericRegrid.h"
//#include "ConvolutionGridder.h"

// Pixel strategies
#include "MedianFilter.h"
#include "MedianFilter_OMP.h"
//#include "NUFFT.h"
//#include "CGSENSE.h"

// Raw & Pixel strategies
#include "DummyRecon.h"
#include "InvertOrder.h"
#include "DumpToFile.h"


ReconContext::ReconContext () {

	m_strategies.begin();

	m_strategies.push_back( ((ReconStrategy*) new DummyRecon())       );
	m_strategies.push_back( ((ReconStrategy*) new HannWindow())       );
	m_strategies.push_back( ((ReconStrategy*) new MedianFilter())     );
	m_strategies.push_back( ((ReconStrategy*) new InvertOrder())      );
	m_strategies.push_back( ((ReconStrategy*) new DumpToFile())       );
    #ifndef __WIN32__
	m_strategies.push_back( ((ReconStrategy*) new MedianFilter_OMP()) );
	#endif
		
}

