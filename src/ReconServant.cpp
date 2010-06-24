#include "ReconServant.h"
#include "ReconContext.h"


/**************************************************************************************************/
ReconServant::ReconServant  ()               {

	m_have_raw    = false;
	m_have_helper = false;
	m_have_pixel  = false;
	m_have_labels = false;
}


/**************************************************************************************************/
ReconServant::~ReconServant ()               {

}


/**************************************************************************************************/
error_code
ReconServant::process_data  (method m)       {

	error_code e = OK;

	ReconContext* context = new ReconContext();
	context->Strategy(m);

	if (m_have_raw)
		context->Strategy()->SetRaw(&m_raw);
	if (m_have_helper)
		context->Strategy()->SetHelper(&m_helper);
	if (m_have_pixel)
		context->Strategy()->SetPixel(&m_pixel);
	if (m_have_labels)
		context->Strategy()->SetLabels(m_labels);

	e = context->ProcessData();
	
	cout << "Finished processing. Getting results ..." << endl;

	if (m_have_raw)
		context->Strategy()->GetRaw(&m_raw);
	if (m_have_helper)
		context->Strategy()->GetHelper(&m_helper);
	if (m_have_pixel)
		context->Strategy()->GetPixel(&m_pixel);
	if (m_have_labels)
		context->Strategy()->GetLabels(m_labels);

	cout << "... done. Will handle control back to client." << endl;
	delete context;

	return e;

	
}


/**************************************************************************************************/
void
ReconServant::raw          (const raw_data& d)   {
	m_raw = d;
	m_have_raw = true;
}


/**************************************************************************************************/
raw_data*
ReconServant::raw          ()                    {
	return new raw_data (m_raw);
}


/**************************************************************************************************/
void
ReconServant::helper       (const raw_data& d)   {
	m_helper = d;
	m_have_helper = true;
}


/**************************************************************************************************/
raw_data*
ReconServant::helper       ()                    {
	return new raw_data (m_helper);
}


/**************************************************************************************************/
void
ReconServant::pixel        (const pixel_data& d) {
	m_pixel = d;
	m_have_pixel = true;
}


/**************************************************************************************************/
pixel_data*
ReconServant::pixel        ()                    {
	return new pixel_data (m_pixel);
}


/**************************************************************************************************/
void 
ReconServant::labels       (const strings& d)    {
	m_labels = d;
	m_have_labels = true;
}


/**************************************************************************************************/
strings* 
ReconServant::labels       ()                    {
		return new strings (m_labels);
}
