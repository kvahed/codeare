
///////////////////////////////////////////////////////////////////////////////
//
//                     Smart Pointers - version 1.9
//
//      For online documentation and newer versions of the code check:
//         http://www.csc.lsu.edu/~kosmas/smartPointers.html
//
//
//   Copyright (C) 1998-2003 Kosmas Karadimitriou, kosmas@computer.org
// 
// You have my permission to use this software in any way you like. To the
// best of my knowledge it works fine, however you use it at your own risk.
// For questions, suggestions, bug reports, feedback, or improvements regarding
// this code please drop me an email - I'd love to hear from you!
//
//----------------------------------------------------------------------------------
// 
//                                DOCUMENTATION
//
// This file contains the implementation of Ptr, an "intrusive" smart pointer that can
// block the most common memory-related problems such as memory leaks, dereferencing
// NULL or unitialized pointers, or trying to access deleted objects. In addition,
// it contains the implementation of AutoPtr, a very useful auto-deletion pointer
// that automatically deletes its contents when it goes out of scope (excellent for
// avoiding all that exception clean up code when allocating temporary memory).
//
// A big bonus of the current implementation is that you can turn the smart pointers off
// by simply commenting out a single line (the _SMART_PTR_ definition). Your code will
// automatically switch to using "almost raw" pointers. Extremely useful for performance
// testing, or for simply turning off the smart pointers in the release version of your
// software... (personally I am against this, since I believe we should always keep the
// "safety net" of smart pointers even in release mode, however I realize that some
// applications may want to squeeze every tiny bit of performance, and going from
// smart pointers back to raw pointers could result in some performance increase).
//
//
// To use smart pointers in your code:
//
// 1. Include this header file
// 2. Replace pointer declarations   CSomeClass*   with   Ptr<CSomeClass>
// 3. For the classes that you want to use smart pointers on (as in the above example
//    CSomeClass), inherit from SmartObject and then define a GetClassName() function,
//    for example:
//
//        class CSomeClass : public SmartObject {
//            ...
//            const char* GetClassName() const { return "CSomeClass"; }
//            ...
//        }
//
//
// When using Ptr, common memory errors would result in a meaningful exception
// instead of crashing your system, or (worse?) go on undetected, e.g.
//
//   Ptr<String> p;              //automatically initialized to NULL
//   p->MakeUpper();             //throws exception (unitialized pointer)
//   p = NEW(String); p = NULL;  //throws exception (memory leak)
//   p = NULL; p->MakeUpper();   //throws exception (NULL pointer)
//   p = NEW(String); delete p; p->MakeUpper(); //throws exception (deleted pointer)
//
// Note that the above errors are detected even if they happen indirectly, e.g.:
//
//   Ptr<String> p1, p2;
//   p1 = NEW(String);
//   p2 = p1;
//   delete p1;
//   p2->MakeUpper();  //throws exception (tried to access deleted object)
//
// If you use "new" instead of "NEW" you won't have memory leak detection,
// but you are still protected against all the other errors mentioned above. 
// Smart pointers Ptr can be interchanged or mixed freely with raw pointers,
// however keep in mind that operations with raw pointers are unprotected.
// The best advise is to always use smart pointers - however, if you do mix
// smart and raw pointers, then at least make sure you delete your objects
// through smart pointers, so that all smart pointers that point to the
// same (deleted) object are updated.
// 
//----------------------------------------------------------------------------------
// 
//                              DISCLAIMER
// 
// This software is provided by the author "as is" and any express or implied
// warranties, including, but not limited to, the implied warranties of 
// merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the author be liable for any direct, indirect, incidental,
// special, exemplary, or consequential damages (including, but not limited to,
// procurement of substitute goods or services; loss of use, data, or profits; 
// or business interruption) however caused and on any theory of liability,
// whether in contract, strict liability, or tort (including negligence or 
// otherwise) arising in any way out of the use of this software.
//
///////////////////////////////////////////////////////////////////////////////


#ifndef _PTR_HEADER_
#define _PTR_HEADER_

#include <iostream>

#undef  _SMART_PTR_

// Commenting out the next line turns smart pointers into "almost raw" pointers,
// when function inlining is enabled (e.g. in Release mode of VC++).
#define _SMART_PTR_  1


// Customize these as needed for your application
#define USER_MSG(msg)  { std::cerr << msg << std::endl; }
#define TRAP(cond, msg) { if(cond) { USER_MSG(msg); throw msg; } }
#define ERROR(msg) { TRAP(true, msg); }


// Note: for the objects that you use NEW instead of new, you get memory leak
// detection. If you don't use NEW is OK, but the memory leak detection for 
// that object will be turned off. Also note that it works only for SmartObjects.
#undef NEW
#define NEW(obj)  New(new obj)

template<class ELT>
ELT New(ELT ptr) 
{
    TRAP(ptr == NULL, "Failed to allocate memory");
    ptr->TrackObjectMemory() = true;
    return ptr;
}



//-----------------------------------------------------------------------------------
// GenericProxy is simply a placeholder for Proxies of different types, since Proxy
// is a template.
//-----------------------------------------------------------------------------------
class GenericProxy {
public:
    GenericProxy() {}
    virtual ~GenericProxy() {}
    virtual void Delete() = 0;
};


//-----------------------------------------------------------------------------------
// Proxy is a helper class to realize the smart pointer scheme. Every SmartObject
// has exactly one proxy. All pointers to a SmartObject are actually pointing
// to its proxy, which remains alive even after the object has been deleted. This
// helps in protecting against trying to access deleted objects, and also is used
// for memory leak detection. The proxy is deleted automatically when no more smart 
// pointers point to it. The user is not actually aware of the proxies at all.
//-----------------------------------------------------------------------------------
template <class ELT> 
class Proxy : public GenericProxy {
public:
    Proxy() { 
        m_ptr = NULL; 
        m_refs = 0; 
    }
    Proxy(ELT* ptr) { 
        m_ptr = ptr; 
        m_refs = 0; 
    }
    virtual ~Proxy() {
        try {
            if (m_ptr != NULL) {
                if (m_ptr->TrackObjectMemory()) {
                    USER_MSG("Memory leak detected! Object type: " << m_ptr->GetClassName());
                }
                m_ptr->SetProxy(NULL);
            }
        }
        catch(...) {
            USER_MSG("Exception in Proxy dtor!");
        }
    }
    operator ELT*() const {
        return m_ptr; 
    }
    ELT& operator*() const {
        TRAP(m_ptr == NULL, "Tried to access deleted object"); 
        return *m_ptr; 
    }
    ELT* operator->() const {
        TRAP(m_ptr == NULL, "Tried to access deleted object"); 
        return m_ptr; 
    }
    void Delete() {
        TRAP(m_ptr == NULL, "Proxy cannot delete NULL pointer"); 
        m_ptr = NULL; 
    }
    void UpCount() { 
        ++m_refs; 
    }
    void DownCount() { 
        if (--m_refs == 0) {
            delete this;
        }
    }
private:
    ELT* m_ptr;
    long m_refs;
    Proxy<ELT>& operator=(Proxy<ELT>& other); //disallow operator=
    Proxy(const Proxy<ELT>& other); //disallow copy ctor
};



#ifdef _SMART_PTR_

//-------------------------------------------------------------------------
// Ptr is the smart pointer class that replaces the "raw" C/C++ pointers.
// For usage comments, see documentation in the beginning of this file.
//-------------------------------------------------------------------------
template <class ELT> 
class Ptr {
public:
    Ptr() : m_targetProxy(NULL) {}
    Ptr(ELT* ptr) {
        if (ptr == NULL) {
            m_targetProxy = NULL;
            return;
        }
        if (ptr->GetProxy() == NULL) {
            //this is the first Ptr who accesses 
            //the object so it must create its proxy
            ptr->SetProxy(new Proxy<ELT>(ptr));
        }
        m_targetProxy = (Proxy<ELT>*)(ptr->GetProxy());
        m_targetProxy->UpCount(); 
    }
    template<class OTHER>
    Ptr(const Ptr<OTHER>& other) {
        ELT* foo = other; //this produces compilation error on illegal pointer conversions
        m_targetProxy = (Proxy<ELT>*)(other.GetTargetProxy());
        if (m_targetProxy != NULL) {
            m_targetProxy->UpCount(); 
        }
    }
    Ptr(const Ptr<ELT>& other) {
        m_targetProxy = other.m_targetProxy;
        if (m_targetProxy != NULL) {
            m_targetProxy->UpCount();
        }
    }
    virtual ~Ptr() { 
        try {
            if (m_targetProxy != NULL) {
                m_targetProxy->DownCount(); 
            }
        }
        catch(...) {
            USER_MSG("Exception in Ptr dtor!");
        }
    }
    operator ELT*() const {
        if (m_targetProxy == NULL) {
            return NULL;
        }
        return m_targetProxy->operator ELT*();
    }
    ELT& operator*() const { 
        TRAP(m_targetProxy == NULL, "Pointer is NULL"); 
        return m_targetProxy->operator*(); 
    }
    ELT* operator->() const { 
        TRAP(m_targetProxy == NULL, "Pointer is NULL"); 
        return m_targetProxy->operator->(); 
    }
    virtual Ptr<ELT>& operator=(Ptr<ELT> &ptr) { 
        return operator=((ELT*)ptr); 
    }
    virtual Ptr<ELT>& operator=(ELT* ptr) {
        if (m_targetProxy != NULL) {
            m_targetProxy->DownCount();
        }
        if (ptr == NULL) {
            m_targetProxy = NULL;
        }
        else {
            if (ptr->GetProxy() == NULL) {
                //this is the first Ptr who accesses 
                //the object so it must create its proxy
                ptr->SetProxy(new Proxy<ELT>(ptr));
            }
            m_targetProxy = (Proxy<ELT>*)(ptr->GetProxy());
            m_targetProxy->UpCount();
        }
        return *this; 
    }
    Proxy<ELT>* GetTargetProxy() const { 
        return m_targetProxy; 
    }
private:
    Proxy<ELT>* m_targetProxy;
};

#else

//-------------------------------------------------------------------------
// This is the fast version of class Ptr, or the "almost raw" pointers.
// In this version Ptr has lost all its smartness, except for being always 
// NULL before initialization, to preserve the semantics of smart pointers.
//-------------------------------------------------------------------------
template <class ELT> 
class Ptr {
public:
    Ptr() : m_ptr(NULL) {}
    Ptr(ELT* ptr) : m_ptr(ptr) {}
    template<class OTHER>
    Ptr(const Ptr<OTHER>& other) {
        m_ptr = other.GetTargetObject();
    }
    virtual ~Ptr() {}
    operator ELT*() const {
        return m_ptr;
    }
    ELT& operator*() const { 
        return *m_ptr;
    }
    ELT* operator->() const { 
        return m_ptr;
    }
    virtual Ptr<ELT>& operator=(Ptr<ELT> &ptr) { 
        m_ptr = ptr;
        return *this; 
    }
    virtual Ptr<ELT>& operator=(ELT* ptr) {
        m_ptr = ptr;
        return *this; 
    }
    ELT* GetTargetObject() const {
        return m_ptr;
    }
private:
    ELT* m_ptr;
};


#endif


//----------------------------------------------------------------------------------
// AutoPtr deletes the memory it points to when it goes out of scope.
// Unlike std::auto_ptr, you CANNOT assign anything to an AutoPtr, you
// can only initialize it to something, e.g: AutoPtr<String> p = NEW(String);
// This is to save us from all the troubles about ownership transfers.
//----------------------------------------------------------------------------------
template <class ELT> 
class AutoPtr : public Ptr<ELT> {
public:
    AutoPtr() : Ptr<ELT>(), m_hasOwnership(false) {}
    AutoPtr(ELT* ptr) : Ptr<ELT>(ptr), m_hasOwnership(true) {}
    AutoPtr(Ptr<ELT>& ptr) : Ptr<ELT>(ptr), m_hasOwnership(true) {}
    virtual ~AutoPtr() { 
        try {
            if (m_hasOwnership) {
                delete *this;
            }
        }
        catch (...) {
            USER_MSG("Exception in AutoPtr dtor!");
        }
    }
private:
    bool m_hasOwnership;
    AutoPtr(const AutoPtr<ELT>& other); //disallow copy ctor
    AutoPtr& operator=(AutoPtr<ELT> &other); //disallow operator=
    virtual Ptr<ELT>& operator=(Ptr<ELT> &ptr) { ERROR("Operator= is disallowed for AutoPtr"); return *this; }
    virtual Ptr<ELT>& operator=(ELT* ptr)  { ERROR("Operator= is disallowed for AutoPtr");  return *this; }
};



//----------------------------------------------------------------------------------
// SmartObject is "smart" in three ways:
//   a) it knows when it is heap allocated (you must use NEW to get this feature)
//   b) it knows and it can return its class name
//   c) it owns the proxy that is being used by smart pointers
// Note that you can use smart pointers only on classes that inherit from SmartObject.
//----------------------------------------------------------------------------------
class SmartObject {
public:
    SmartObject() { 
        m_proxy = NULL;
        m_trackObjectMemory = false;
    }
    virtual ~SmartObject() { 
        try {
            if (m_proxy != NULL) {
                m_proxy->Delete(); 
                m_proxy = NULL; 
            }
        }
        catch (...) {
            USER_MSG("Exception in SmartObject dtor!");
        }
    }
    SmartObject& operator=(const SmartObject& other) { 
        //don't copy anything - every object should have its own proxy
        return *this; 
    }
    //CAUTION: when defining a copy ctor in a SmartObject child, don't forget to add:
    // : SmartObject(other)
    SmartObject(const SmartObject& other) { 
        m_proxy = NULL;
        m_trackObjectMemory = false;
    }
    GenericProxy* GetProxy() const {
        return m_proxy;
    }
    void SetProxy(GenericProxy* proxy) {
        TRAP(m_proxy != NULL && proxy != NULL, "Object has already a proxy");
        m_proxy = proxy;
    }
    bool& TrackObjectMemory() {
        return m_trackObjectMemory;
    }
    virtual const char* GetClassName() const = 0;
private:
    GenericProxy* m_proxy;
    bool m_trackObjectMemory;
};


#endif

