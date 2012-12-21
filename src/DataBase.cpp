#include "DataBase.hpp"

#if defined(__APPLE__)
#  include "AppleDigest.hpp"
#else
#  include "Digest.hpp"
#endif

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
	
  	while (!m_cxfl.empty()) {
		delete m_cxfl.begin()->second;
		m_cxfl.erase(m_cxfl.begin());
	}
			   
  	while (!m_cxdb.empty()) {
		delete m_cxdb.begin()->second;
		m_cxdb.erase(m_cxdb.begin());
	}
			   
  	while (!m_rlfl.empty()) {
		delete m_rlfl.begin()->second;
		m_rlfl.erase(m_rlfl.begin());
	}
	
  	while (!m_rldb.empty()) {
		delete m_rldb.begin()->second;
		m_rldb.erase(m_rldb.begin());
	}
	
  	while (!m_shrt.empty()) {
		delete m_shrt.begin()->second;
		m_shrt.erase(m_shrt.begin());
	}
	
  	while (!m_long.empty()) {
		delete m_long.begin()->second;
		m_long.erase(m_long.begin());
	}
	
  	while (!m_ref.empty())
		m_ref.erase(m_ref.begin());

	return OK;
	
}


template<> Matrix<cxfl>& 
DataBase::AddMatrix (const string name, Ptr< Matrix<cxfl> > m) {
	
	std::string tag[2];
	tag[1] = sha256(name);
	tag[0] = "cxfl";

	assert (m_ref.find (name) == m_ref.end());
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_cxfl.insert (pair< string, Ptr < Matrix<cxfl> > >(tag[0], m));

	return *m_cxfl[tag[0]];
	
}

	
template<> Matrix<cxfl>& 
DataBase::Get<cxfl> (const string name) {
	
	map<string,string[2]>::iterator it = m_ref.find(name);

	return *m_cxfl[it->second[0]];
	
}


template<> void
DataBase::GetMatrix (const string name, cxfl_data& c) {
	
	if (m_ref.find (name) == m_ref.end())
		return;
		
	map<string,string[2]>::iterator it = m_ref.find(name);

	Ptr< Matrix<cxfl> > tmp = m_cxfl[it->second[0]];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		c.dims[j] = tmp->Dim(j);
		c.res[j]  = tmp->Res(j);
	}
	
	c.vals.length(2 * tmp->Size()); 
	
	memcpy (&c.vals[0], &tmp[0], c.vals.length() * sizeof(double));
	
}	


template<> void
DataBase::GetMatrix (const string name, Matrix<cxfl>& m) {
	
	if (m_ref.find (name) == m_ref.end())
		return;

	map<string,string[2]>::iterator it = m_ref.find(name);
		
	m = *m_cxfl[it->second[0]];
	
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
	
	if (m_ref.find (name) != m_ref.end())
		Free<cxfl> (name);

	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "cxfl";
	
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_cxfl.insert (pair<string, Ptr< Matrix<cxfl> > > (tag[0], tmp = NEW (Matrix<cxfl>(mdims, mress))));
	
	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp[0], &c.vals[0], c.vals.length() * sizeof(float));
	
}


template<> void 
DataBase::SetMatrix (const string name, Matrix<cxfl>& m)   {
	
	string tag[2];
	tag[0] = sha256(name);
	tag[1] = "cxfl";

	if (m_ref.find (name) == m_ref.end()) {

		m_ref.insert (pair<string,string[2]> (name,tag));
		m_cxfl.insert (pair<string, Ptr< Matrix<cxfl> > > (tag[0], NEW (Matrix<cxfl>())));

	}
	
	m.SetClassName(name.c_str());

	*m_cxfl[tag[0]] = m;
	
}


template<> Matrix<cxdb>& 
DataBase::AddMatrix (const string name, Ptr< Matrix<cxdb> > m) {
	
	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "cxdb";

	assert (m_ref.find (name) == m_ref.end());
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_cxdb.insert (pair< string, Ptr < Matrix<cxdb> > >(tag[0], m));

	return *m_cxdb[tag[0]];

}

	
template<> Matrix<cxdb>& 
DataBase::Get<cxdb> (const string name) {
	
	map<string,string[2]>::iterator it = m_ref.find(name);

	return *m_cxdb[it->second[0]];
	
}
	

template<> void
DataBase::GetMatrix (const string name, cxdb_data& c) {
	
	if (m_ref.find (name) == m_ref.end())
		return;
		
	map<string,string[2]>::iterator it = m_ref.find(name);

	Ptr< Matrix<cxdb> > tmp = m_cxdb[it->second[0]];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		c.dims[j] = tmp->Dim(j);
		c.res[j]  = tmp->Res(j);
	}
	
	c.vals.length(2 * tmp->Size()); 
	
	memcpy (&c.vals[0], &tmp[0], c.vals.length() * sizeof(double));
	
}	


template<> void
DataBase::GetMatrix (const string name, Matrix<cxdb>& m) {
	
	if (m_ref.find (name) == m_ref.end())
		return;

	map<string,string[2]>::iterator it = m_ref.find(name);
		
	m = *m_cxdb[it->second[0]];
	
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
	
	if (m_ref.find (name) != m_ref.end())
		Free<cxdb> (name);

	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "cxdb";
	
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_cxdb.insert (pair<string, Ptr< Matrix<cxdb> > > (tag[0], tmp = NEW (Matrix<cxdb>(mdims, mress))));

	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp[0], &c.vals[0], c.vals.length() * sizeof(double));
	
}


template<> void 
DataBase::SetMatrix (const string name, Matrix<cxdb>& m)   {
	

	string tag[2];
	tag[0] = sha256(name);
	tag[1] = "cxdb";
	
	if (m_ref.find (name) == m_ref.end()) {

		m_ref.insert (pair<string,string[2]> (name,tag));
		m_cxdb.insert (pair<string, Ptr< Matrix<cxdb> > > (tag[0], NEW (Matrix<cxdb>())));

	}
	
	m.SetClassName(name.c_str());
	
	*m_cxdb[tag[0]] = m;
	
}


template<> Matrix<float>&
DataBase::AddMatrix         (const string name, Ptr< Matrix<float> > m) {
	
	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "float";

	assert (m_ref.find (name) == m_ref.end());
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_rlfl.insert (pair< string, Ptr < Matrix<float> > >(tag[0], m));

	return *m_rlfl[tag[0]];

}
		

template<> Matrix<float>& 
DataBase::Get<float> (const string name) {
	
	map<string,string[2]>::iterator it = m_ref.find(name);

	return *m_rlfl[it->second[0]];
	
}

		
template<> void 
DataBase::GetMatrix        (const string name, rlfl_data& r) {
	
	if (m_ref.find (name) == m_ref.end())
		return;

	map<string,string[2]>::iterator it = m_ref.find(name);

	Ptr< Matrix<float> > tmp = m_rlfl[it->second[0]];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		r.dims[j] = tmp->Dim(j);
		r.res[j]  = tmp->Res(j);
	}
	
	r.vals.length(tmp->Size()); 
	
	memcpy (&r.vals[0], &tmp[0], tmp->Size() * sizeof(float));
	
}
		

template<> void
DataBase::GetMatrix         (const string name, Matrix<float>& m) {
	
	if (m_ref.find (name) == m_ref.end())
		return;
		
	map<string,string[2]>::iterator it = m_ref.find(name);

	m = *m_rlfl[it->second[0]];
	
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
	
	if (m_ref.find (name) != m_ref.end())
		Free<float> (name);

	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "float";
	
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_rlfl.insert (pair<string, Ptr< Matrix<float> > > (tag[0], tmp = NEW (Matrix<float>(mdims, mress))));

	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp[0], &r.vals[0], tmp->Size() * sizeof(float));
	
}
		
		
template<> void 
DataBase::SetMatrix         (const string name, Matrix<float>& m)   {
	
	string tag[2];
	tag[0] = sha256(name);
	tag[1] = "float";

	if (m_ref.find (name) == m_ref.end()) {

		m_ref.insert (pair<string,string[2]> (name,tag));
		m_rlfl.insert (pair<string, Ptr< Matrix<float> > > (tag[0], NEW (Matrix<float>())));

	}
	
	m.SetClassName(name.c_str());
	
	*m_rlfl[tag[0]] = m;
	
}
		

template<> Matrix<double>&
DataBase::AddMatrix         (const string name, Ptr< Matrix<double> > m) {
	
	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "double";

	assert (m_ref.find (name) == m_ref.end());
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_rldb.insert (pair< string, Ptr < Matrix<double> > >(tag[0], m));

	return *m_rldb[tag[0]];

}
		

template<> Matrix<double>& 
DataBase::Get<double> (const string name) {
			
	map<string,string[2]>::iterator it = m_ref.find(name);

	return *m_rldb[it->second[0]];
	
}

		
template<> void 
DataBase::GetMatrix        (const string name, rldb_data& r) {
	
	if (m_ref.find (name) == m_ref.end())
		return;
		
	map<string,string[2]>::iterator it = m_ref.find(name);

	Ptr< Matrix<double> > tmp = m_rldb[it->second[0]];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		r.dims[j] = tmp->Dim(j);
		r.res[j]  = tmp->Res(j);
	}
	
	r.vals.length(tmp->Size()); 
	
	memcpy (&r.vals[0], &tmp[0], tmp->Size() * sizeof(double));
	
}
		

template<> void
DataBase::GetMatrix         (const string name, Matrix<double>& m) {
	
	if (m_ref.find (name) == m_ref.end())
		return;

	map<string,string[2]>::iterator it = m_ref.find(name);

	m = *m_rldb[it->second[0]];
	
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
	
	if (m_ref.find (name) != m_ref.end())
		Free<double> (name);

	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "double";
	
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_rldb.insert (pair<string, Ptr< Matrix<double> > > (tag[0], tmp = NEW (Matrix<double>(mdims, mress))));

	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp[0], &r.vals[0], tmp->Size() * sizeof(double));
	
}
		
		
template<> void 
DataBase::SetMatrix         (const string name, Matrix<double>& m)   {
	
	string tag[2];
	tag[0] = sha256(name);
	tag[1] = "double";
	
	if (m_ref.find (name) == m_ref.end()) {

		m_ref.insert (pair<string,string[2]> (name,tag));
		m_rldb.insert (pair<string, Ptr< Matrix<double> > > (tag[0], NEW (Matrix<double>())));

	}
	
	m.SetClassName(name.c_str());
	
	*m_rldb[tag[0]] = m;	
	
}
	

template<> Matrix<short>&
DataBase::AddMatrix        (const string name, Ptr< Matrix<short> > m) {
	
	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "short";

	assert (m_ref.find (name) == m_ref.end());
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_shrt.insert (pair< string, Ptr < Matrix<short> > >(tag[0], m));

	return *m_shrt[tag[0]];

}


template<> Matrix<short>&   
DataBase::Get<short>         (const string name) {
	
	map<string,string[2]>::iterator it = m_ref.find(name);

	return *m_shrt[it->second[0]];
	
}


template<> void 
DataBase::GetMatrix          (const string name, shrt_data& p) {
	
	if (m_ref.find (name) == m_ref.end())
		return;
		
	map<string,string[2]>::iterator it = m_ref.find(name);

	Ptr< Matrix<short> > tmp = m_shrt[it->second[0]];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		p.dims[j] = tmp->Dim(j);
		p.res[j]  = tmp->Res(j);
	}
	
	p.vals.length(tmp->Size()); 
	
	memcpy (&p.vals[0], &tmp[0], tmp->Size() * sizeof(short));
	
	Free<short>(name);
	
}


template<> void
DataBase::GetMatrix         (const string name, Matrix<short>& m) {
		
	if (m_ref.find (name) == m_ref.end())
		return;
		
	map<string,string[2]>::iterator it = m_ref.find(name);

	m = *m_shrt[it->second[0]];
	
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
	
	if (m_ref.find (name) != m_ref.end())
		Free<short> (name);

	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "short";
	
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_shrt.insert (pair<string, Ptr< Matrix<short> > > (tag[0], tmp = NEW (Matrix<short>(mdims, mress))));

	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp[0], &p.vals[0], tmp->Size() * sizeof(short));
	
}


template<> void 
DataBase::SetMatrix         (const string name, Matrix<short>& m)   {
	
	string tag[2];
	tag[0] = sha256(name);
	tag[1] = "short";
		
	if (m_ref.find (name) == m_ref.end()) {

		m_ref.insert (pair<string,string[2]> (name,tag));
		m_shrt.insert (pair<string, Ptr< Matrix<short> > > (tag[0], NEW (Matrix<short>())));

	}
	
	m.SetClassName(name.c_str());
	
	*m_shrt[tag[0]] = m;	

}

		
template<> Matrix<long>&
DataBase::AddMatrix        (const string name, Ptr< Matrix<long> > m) {
	
	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "long";

	assert (m_ref.find (name) == m_ref.end());
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_long.insert (pair< string, Ptr < Matrix<long> > >(tag[0], m));

	return *m_long[tag[0]];

}


template<> Matrix<long>&   
DataBase::Get<long>         (const string name) {
	
	map<string,string[2]>::iterator it = m_ref.find(name);

	return *m_long[it->second[0]];
	
}


template<> void 
DataBase::GetMatrix          (const string name, long_data& p) {
	
	if (m_ref.find (name) == m_ref.end())
		return;
		
	map<string,string[2]>::iterator it = m_ref.find(name);

	Ptr< Matrix<long> > tmp = m_long[it->second[0]];
	
	for (int j = 0; j < INVALID_DIM; j++) {
		p.dims[j] = tmp->Dim(j);
		p.res[j]  = tmp->Res(j);
	}
	
	p.vals.length(tmp->Size()); 
	
	memcpy (&p.vals[0], &tmp[0], tmp->Size() * sizeof(long));
	
	Free<long> (name);
	
}


template<> void
DataBase::GetMatrix         (const string name, Matrix<long>& m) {
		
	if (m_ref.find (name) == m_ref.end())
		return;
		
	map<string,string[2]>::iterator it = m_ref.find(name);

	m = *m_long[it->second[0]];
	
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
	
	if (m_ref.find (name) != m_ref.end())
		Free<long> (name);

	std::string tag[2];
	tag[0] = sha256(name);
	tag[1] = "long";
	
	m_ref.insert (pair<string, string[2]> (name, tag));
	m_long.insert (pair<string, Ptr< Matrix<long> > > (tag[0], tmp = NEW (Matrix<long>(mdims, mress))));

	tmp->SetClassName(name.c_str());
	
	memcpy (&tmp[0], &p.vals[0], tmp->Size() * sizeof(long));
	
}


template<> void 
DataBase::SetMatrix         (const string name, Matrix<long>& m)   {
	
	string tag[2];
	tag[0] = sha256(name);
	tag[1] = "long";
		
	if (m_ref.find (name) == m_ref.end()) {

		m_ref.insert (pair<string,string[2]> (name,tag));
		m_long.insert (pair<string, Ptr< Matrix<long> > > (tag[0], NEW (Matrix<long>())));

	}
	
	m.SetClassName(name.c_str());
	
	*m_long[tag[0]] = m;	

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
