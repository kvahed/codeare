#include "DataBase.hpp"

DataBase* DataBase::m_inst = 0; 

DataBase::~DataBase () { 
	
	Finalise();
	m_inst = 0; 
	
}


DataBase* 
DataBase::Instance ()  {

    if (m_inst == 0) 
        m_inst = new DataBase ();

    return m_inst;
		
}
	
error_code
DataBase::Finalise () {
	
	while (!m_cxfl.empty())
		Free<cxfl>(m_cxfl.begin()->first.c_str());
	
	while (!m_cxdb.empty()) 
		Free<cxdb>(m_cxdb.begin()->first.c_str());
	
	while (!m_rlfl.empty())
		Free<float>(m_rlfl.begin()->first.c_str());
	
	while (!m_rldb.empty()) 
		Free<double>(m_rldb.begin()->first.c_str());
	
	while (!m_shrt.empty())
		Free<short>(m_shrt.begin()->first.c_str());
	
	while (!m_long.empty())
		Free<long>(m_long.begin()->first.c_str());
	
	return OK;
	
}


template<> Matrix<cxfl>& 
DataBase::AddMatrix (const string name, Ptr< Matrix<cxfl> > m) {
	
	assert (m_cxfl.find (name) == m_cxfl.end());
	m_cxfl.insert (pair< string, Ptr < Matrix<cxfl> > >(name, m));
	return *m_cxfl[name];
	
}

	
template<> Matrix<cxfl>& 
DataBase::Get<cxfl> (const string name) {
	
	return *m_cxfl[name];
	
}
	

template<> void
DataBase::GetMatrix (const string name, cxfl_data& c) {

	if (m_cxfl.find (name) == m_cxfl.end())
		return;
		
	Ptr< Matrix<cxfl> > tmp = m_cxfl[name];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		c.dims[j] = tmp->Dim(j);
		c.res[j]  = tmp->Res(j);
	}
	
	c.vals.length(2 * tmp->Size()); 
	
	memcpy (&c.vals[0], &tmp->At(0), c.vals.length() * sizeof(float));
	
}	


template<> void
DataBase::GetMatrix (const string name, Matrix<cxfl>& m) {
	
	if (m_cxfl.find (name) == m_cxfl.end())
		return;
	
	m = *m_cxfl[name];
	
}


template<> void 
DataBase::SetMatrix (const string name, const cxfl_data& c)   {

	size_t mdims [INVALID_DIM];
	float  mress [INVALID_DIM];
	
	for (int i = 0; i < INVALID_DIM; i++) {
		mdims[i] = c.dims[i];
		mress[i] = c.res[i];
	}
	
	Ptr< Matrix<cxfl> > tmp;
	
	if (m_cxfl.find (name) != m_cxfl.end())
		Free<cxfl> (name);

	m_cxfl.insert (pair<string, Ptr< Matrix<cxfl> > > (name, tmp = NEW (Matrix<cxfl>(mdims, mress))));
	
	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp->At(0), &c.vals[0], c.vals.length() * sizeof(float));
	
}


template<> void 
DataBase::SetMatrix (const string name, Matrix<cxfl>& m)   {
	
	if (m_cxfl.find (name) == m_cxfl.end())
		m_cxfl.insert (pair<string, Ptr< Matrix<cxfl> > > (name, NEW (Matrix<cxfl>())));
	
	m.SetClassName(name.c_str());
	
	*m_cxfl[name] = m;
	
}


template<> Matrix<cxdb>& 
DataBase::AddMatrix (const string name, Ptr< Matrix<cxdb> > m) {
	
	assert (m_cxdb.find (name) == m_cxdb.end());
	m_cxdb.insert (pair< string, Ptr < Matrix<cxdb> > >(name, m));
	return *m_cxdb[name];
	
}

	
template<> Matrix<cxdb>& 
DataBase::Get<cxdb> (const string name) {
	
	return *m_cxdb[name];
	
}
	

template<> void
DataBase::GetMatrix (const string name, cxdb_data& c) {
	
	if (m_cxdb.find (name) == m_cxdb.end())
		return;
		
	Ptr< Matrix<cxdb> > tmp = m_cxdb[name];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		c.dims[j] = tmp->Dim(j);
		c.res[j]  = tmp->Res(j);
	}
	
	c.vals.length(2 * tmp->Size()); 
	
	memcpy (&c.vals[0], &tmp->At(0), c.vals.length() * sizeof(double));
	
}	


template<> void
DataBase::GetMatrix (const string name, Matrix<cxdb>& m) {
	
	if (m_cxdb.find (name) == m_cxdb.end())
		return;
	
	m = *m_cxdb[name];
	
}


template<> void 
DataBase::SetMatrix (const string name, const cxdb_data& c)   {
	
	size_t mdims [INVALID_DIM];
	float  mress [INVALID_DIM];
	
	for (int i = 0; i < INVALID_DIM; i++) {
		mdims[i] = c.dims[i];
		mress[i] = c.res[i];
	}
	
	Ptr< Matrix<cxdb> > tmp;
	
	if (m_cxdb.find (name) != m_cxdb.end())
		Free<cxdb> (name);
		
	m_cxdb.insert (pair<string, Ptr< Matrix<cxdb> > > (name, tmp = NEW (Matrix<cxdb>(mdims, mress))));

	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp->At(0), &c.vals[0], c.vals.length() * sizeof(double));
	
}


template<> void 
DataBase::SetMatrix (const string name, Matrix<cxdb>& m)   {
	
	if (m_cxdb.find (name) == m_cxdb.end())
		m_cxdb.insert (pair<string, Ptr< Matrix<cxdb> > > (name, NEW (Matrix<cxdb>())));
	
	m.SetClassName(name.c_str());
	
	*m_cxdb[name] = m;
	
}


template<> Matrix<float>&
DataBase::AddMatrix         (const string name, Ptr< Matrix<float> > m) {
	
	assert(m_rlfl.find (name) == m_rlfl.end());
	m_rlfl.insert (pair< string, Ptr < Matrix<float> > >(name, m));
	return *m_rlfl[name];
	
}
		

template<> Matrix<float>& 
DataBase::Get<float> (const string name) {
			
	return *m_rlfl[name];
	
}

		
template<> void 
DataBase::GetMatrix        (const string name, rlfl_data& r) {
	
	if (m_rlfl.find (name) == m_rlfl.end())
		return;
	
	Ptr< Matrix<float> > tmp = m_rlfl[name];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		r.dims[j] = tmp->Dim(j);
		r.res[j]  = tmp->Res(j);
	}
	
	r.vals.length(tmp->Size()); 
	
	memcpy (&r.vals[0], &tmp->At(0), tmp->Size() * sizeof(float));
	
}
		

template<> void
DataBase::GetMatrix         (const string name, Matrix<float>& m) {
	
	if (m_rlfl.find (name) == m_rlfl.end())
		return;
	
	m = *m_rlfl[name];
	
}
		
	
template<> void
DataBase::SetMatrix        (const string name, const rlfl_data& r)   {
	
	size_t mdims [INVALID_DIM];
	float  mress [INVALID_DIM];
	
	for (int i = 0; i < INVALID_DIM; i++) {
		mdims[i] = r.dims[i];
		mress[i] = r.res[i];
	}
	
	Ptr< Matrix<float> > tmp;
	
	if (m_rlfl.find (name) != m_rlfl.end())
		Free<float> (name);

	m_rlfl.insert (pair<string, Ptr< Matrix<float> > > (name, tmp = NEW( Matrix<float>(mdims, mress))));
	
	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp->At(0), &r.vals[0], tmp->Size() * sizeof(float));
	
}
		
		
template<> void 
DataBase::SetMatrix         (const string name, Matrix<float>& m)   {
	
	if (m_rlfl.find (name) == m_rlfl.end())
		m_rlfl.insert (pair<string, Ptr< Matrix<float> > > (name, NEW( Matrix<float>())));
	
	m.SetClassName(name.c_str());
	
	*m_rlfl[name] = m;
	
}
		

template<> Matrix<double>&
DataBase::AddMatrix         (const string name, Ptr< Matrix<double> > m) {
	
	assert(m_rldb.find (name) == m_rldb.end());
	m_rldb.insert (pair< string, Ptr < Matrix<double> > >(name, m));
	return *m_rldb[name];
	
}
		

template<> Matrix<double>& 
DataBase::Get<double> (const string name) {
			
	return *m_rldb[name];
	
}

		
template<> void 
DataBase::GetMatrix        (const string name, rldb_data& r) {
	
	if (m_rldb.find (name) == m_rldb.end())
		return;
	
	Ptr< Matrix<double> > tmp = m_rldb[name];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		r.dims[j] = tmp->Dim(j);
		r.res[j]  = tmp->Res(j);
	}
	
	r.vals.length(tmp->Size()); 
	
	memcpy (&r.vals[0], &tmp->At(0), tmp->Size() * sizeof(double));
	
}
		

template<> void
DataBase::GetMatrix         (const string name, Matrix<double>& m) {
	
	if (m_rldb.find (name) == m_rldb.end())
		return;
	
	m = *m_rldb[name];
	
}
		
	
template<> void 
DataBase::SetMatrix        (const string name, const rldb_data& r)   {
	
	size_t mdims [INVALID_DIM];
	float  mress [INVALID_DIM];
	
	for (int i = 0; i < INVALID_DIM; i++) {
		mdims[i] = r.dims[i];
		mress[i] = r.res[i];
	}
	
	Ptr< Matrix<double> > tmp;
	
	if (m_rldb.find (name) != m_rldb.end())
		Free<double> (name);

	m_rldb.insert (pair<string, Ptr< Matrix<double> > > (name, tmp = NEW( Matrix<double>(mdims, mress))));
	
	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp->At(0), &r.vals[0], tmp->Size() * sizeof(double));
	
}
		
		
template<> void 
DataBase::SetMatrix         (const string name, Matrix<double>& m)   {
	
	if (m_rldb.find (name) == m_rldb.end())
		m_rldb.insert (pair<string, Ptr< Matrix<double> > > (name, NEW( Matrix<double>())));
	
	m.SetClassName(name.c_str());
	
	*m_rldb[name] = m;
	
}
	

template<> Matrix<short>&
DataBase::AddMatrix        (const string name, Ptr< Matrix<short> > m) {
	
	assert (m_shrt.find (name) == m_shrt.end());
	m_shrt.insert (pair<string, Ptr< Matrix<short> > > (name, m));
	return *m_shrt[name];
	
}


template<> Matrix<short>&   
DataBase::Get<short>         (const string name) {
	
	return *m_shrt[name];
	
}


template<> void 
DataBase::GetMatrix          (const string name, shrt_data& p) {
	
	if (m_shrt.find (name) == m_shrt.end())
		return;
	
	Ptr< Matrix<short> > tmp = m_shrt[name];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		p.dims[j] = tmp->Dim(j);
		p.res[j]  = tmp->Res(j);
	}
	
	p.vals.length(tmp->Size()); 
	
	memcpy (&p.vals[0], &tmp->At(0), tmp->Size() * sizeof(short));
	
	Free<short>(name);
	
}


template<> void
DataBase::GetMatrix         (const string name, Matrix<short>& m) {
		
	if (m_shrt.find (name) == m_shrt.end())
		return;
	
	m = *m_shrt[name];
	
}

		
template<> void 
DataBase::SetMatrix         (const string name, const shrt_data& p)   {
	
	size_t mdims [INVALID_DIM];
	float  mress [INVALID_DIM];
	
	for (int i = 0; i < INVALID_DIM; i++) {
		mdims[i] = p.dims[i];
		mress[i] = p.res[i];
	}
	
	Ptr< Matrix<short> > tmp;
	
	if (m_shrt.find (name) != m_shrt.end())
		Free<short> (name);

	m_shrt.insert (pair<string, Ptr< Matrix<short> > > (name, tmp = NEW (Matrix<short>(mdims, mress))));

	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp->At(0), &p.vals[0], tmp->Size() * sizeof(short));
	
}


template<> void 
DataBase::SetMatrix         (const string name, Matrix<short>& m)   {
	
	if (m_shrt.find (name) == m_shrt.end())
		m_shrt.insert (pair<string, Ptr< Matrix<short> > > (name, NEW (Matrix<short>())));
	
	m.SetClassName(name.c_str());
	
	*m_shrt[name] = m;
	
}

		
template<> Matrix<long>&
DataBase::AddMatrix        (const string name, Ptr< Matrix<long> > m) {
	
	assert (m_long.find (name) == m_long.end());
	m_long.insert (pair<string, Ptr< Matrix<long> > > (name, m));
	return *m_long[name];
	
}


template<> Matrix<long>&   
DataBase::Get<long>         (const string name) {
	
	return *m_long[name];
	
}


template<> void 
DataBase::GetMatrix          (const string name, long_data& p) {
	
	if (m_long.find (name) == m_long.end())
		return;
	
	Ptr< Matrix<long> > tmp = m_long[name];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		p.dims[j] = tmp->Dim(j);
		p.res[j]  = tmp->Res(j);
	}
	
	p.vals.length(tmp->Size()); 
	
	memcpy (&p.vals[0], &tmp->At(0), tmp->Size() * sizeof(long));
	
	Free<long> (name);
	
}


template<> void
DataBase::GetMatrix         (const string name, Matrix<long>& m) {
		
	if (m_long.find (name) == m_long.end())
		return;
	
	m = *m_long[name];
	
}

		
template<> void 
DataBase::SetMatrix         (const string name, const long_data& p)   {
	
	size_t mdims [INVALID_DIM];
	float  mress [INVALID_DIM];
	
	for (int i = 0; i < INVALID_DIM; i++) {
		mdims[i] = p.dims[i];
		mress[i] = p.res[i];
	}
	
	Ptr< Matrix<long> > tmp;
	
	if (m_long.find (name) != m_long.end())
		Free<long> (name);

	m_long.insert (pair<string, Ptr< Matrix<long> > > (name, tmp = NEW (Matrix<long>(mdims, mress))));
	
	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp->At(0), &p.vals[0], tmp->Size() * sizeof(long));
	
}


template<> void 
DataBase::SetMatrix         (const string name, Matrix<long>& m)   {
	
	if (m_long.find (name) == m_long.end())
		m_long.insert (pair<string, Ptr< Matrix<long> > > (name, NEW (Matrix<long>())));
	
	m.SetClassName(name.c_str());
	
	*m_long[name] = m;
	
}

		
map < string, Ptr< Matrix<cxfl> > >& 
DataBase::CXFLMap () {
	return m_cxfl;
}


map < string, Ptr< Matrix<cxdb> > >& 
DataBase::CXDBMap () {
	return m_cxdb;
}


map < string, Ptr< Matrix<float> > >& 
DataBase::RLFLMap () {
	return m_rlfl;
}


map < string, Ptr< Matrix<double> > >& 
DataBase::RLDBMap () {
	return m_rldb;
}


map < string, Ptr< Matrix<short> > >& 
DataBase::SHRTMap () {
	return m_shrt;
}


map < string, Ptr< Matrix<long> > >& 
DataBase::LONGMap () {
	return m_long;
}
