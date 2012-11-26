// Copyright (C) 2005-2010 Code Synthesis Tools CC
//
// This program was generated by CodeSynthesis XSD, an XML Schema to
// C++ data binding compiler.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// In addition, as a special exception, Code Synthesis Tools CC gives
// permission to link this program with the Xerces-C++ library (or with
// modified versions of Xerces-C++ that use the same license as Xerces-C++),
// and distribute linked combinations including the two. You must obey
// the GNU General Public License version 2 in all respects for all of
// the code used other than Xerces-C++. If you modify this copy of the
// program, you may extend this exception to your version of the program,
// but you are not obligated to do so. If you do not wish to do so, delete
// this exception statement from your version.
//
// Furthermore, Code Synthesis Tools CC makes a special exception for
// the Free/Libre and Open Source Software (FLOSS) which is described
// in the accompanying FLOSSE file.
//

#ifndef INPUT_PARAMETERS_HXX
#define INPUT_PARAMETERS_HXX

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/config.hxx>

#if (XSD_INT_VERSION != 3030000L)
#error XSD runtime version mismatch
#endif

#include <xsd/cxx/pre.hxx>

#ifndef XSD_USE_CHAR
#define XSD_USE_CHAR
#endif

#ifndef XSD_CXX_TREE_USE_CHAR
#define XSD_CXX_TREE_USE_CHAR
#endif

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/types.hxx>

#include <xsd/cxx/xml/error-handler.hxx>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

#include <xsd/cxx/tree/parsing.hxx>
#include <xsd/cxx/tree/parsing/byte.hxx>
#include <xsd/cxx/tree/parsing/unsigned-byte.hxx>
#include <xsd/cxx/tree/parsing/short.hxx>
#include <xsd/cxx/tree/parsing/unsigned-short.hxx>
#include <xsd/cxx/tree/parsing/int.hxx>
#include <xsd/cxx/tree/parsing/unsigned-int.hxx>
#include <xsd/cxx/tree/parsing/long.hxx>
#include <xsd/cxx/tree/parsing/unsigned-long.hxx>
#include <xsd/cxx/tree/parsing/boolean.hxx>
#include <xsd/cxx/tree/parsing/float.hxx>
#include <xsd/cxx/tree/parsing/double.hxx>
#include <xsd/cxx/tree/parsing/decimal.hxx>

namespace xml_schema
{
  // anyType and anySimpleType.
  //
  typedef ::xsd::cxx::tree::type type;
  typedef ::xsd::cxx::tree::simple_type< type > simple_type;
  typedef ::xsd::cxx::tree::type container;

  // 8-bit
  //
  typedef signed char byte;
  typedef unsigned char unsigned_byte;

  // 16-bit
  //
  typedef short short_;
  typedef unsigned short unsigned_short;

  // 32-bit
  //
  typedef int int_;
  typedef unsigned int unsigned_int;

  // 64-bit
  //
  typedef long long long_;
  typedef unsigned long long unsigned_long;

  // Supposed to be arbitrary-length integral types.
  //
  typedef long long integer;
  typedef long long non_positive_integer;
  typedef unsigned long long non_negative_integer;
  typedef unsigned long long positive_integer;
  typedef long long negative_integer;

  // Boolean.
  //
  typedef bool boolean;

  // Floating-point types.
  //
  typedef float float_;
  typedef double double_;
  typedef double decimal;

  // String types.
  //
  typedef ::xsd::cxx::tree::string< char, simple_type > string;
  typedef ::xsd::cxx::tree::normalized_string< char, string > normalized_string;
  typedef ::xsd::cxx::tree::token< char, normalized_string > token;
  typedef ::xsd::cxx::tree::name< char, token > name;
  typedef ::xsd::cxx::tree::nmtoken< char, token > nmtoken;
  typedef ::xsd::cxx::tree::nmtokens< char, simple_type, nmtoken > nmtokens;
  typedef ::xsd::cxx::tree::ncname< char, name > ncname;
  typedef ::xsd::cxx::tree::language< char, token > language;

  // ID/IDREF.
  //
  typedef ::xsd::cxx::tree::id< char, ncname > id;
//  typedef ::xsd::cxx::tree::idref< char, ncname, type > idref;
  typedef ::xsd::cxx::tree::idrefs< char, simple_type, idref > idrefs;

  // URI.
  //
  typedef ::xsd::cxx::tree::uri< char, simple_type > uri;

  // Qualified name.
  //
  typedef ::xsd::cxx::tree::qname< char, simple_type, uri, ncname > qname;

  // Binary.
  //
  typedef ::xsd::cxx::tree::buffer< char > buffer;
  typedef ::xsd::cxx::tree::base64_binary< char, simple_type > base64_binary;
  typedef ::xsd::cxx::tree::hex_binary< char, simple_type > hex_binary;

  // Date/time.
  //
  typedef ::xsd::cxx::tree::time_zone time_zone;
  typedef ::xsd::cxx::tree::date< char, simple_type > date;
  typedef ::xsd::cxx::tree::date_time< char, simple_type > date_time;
  typedef ::xsd::cxx::tree::duration< char, simple_type > duration;
  typedef ::xsd::cxx::tree::gday< char, simple_type > gday;
  typedef ::xsd::cxx::tree::gmonth< char, simple_type > gmonth;
  typedef ::xsd::cxx::tree::gmonth_day< char, simple_type > gmonth_day;
  typedef ::xsd::cxx::tree::gyear< char, simple_type > gyear;
  typedef ::xsd::cxx::tree::gyear_month< char, simple_type > gyear_month;
  typedef ::xsd::cxx::tree::time< char, simple_type > time;

  // Entity.
  //
  typedef ::xsd::cxx::tree::entity< char, ncname > entity;
  typedef ::xsd::cxx::tree::entities< char, simple_type, entity > entities;

  // Flags and properties.
  //
  typedef ::xsd::cxx::tree::flags flags;
  typedef ::xsd::cxx::tree::properties< char > properties;

  // Parsing/serialization diagnostics.
  //
  typedef ::xsd::cxx::tree::severity severity;
  typedef ::xsd::cxx::tree::error< char > error;
  typedef ::xsd::cxx::tree::diagnostics< char > diagnostics;

  // Exceptions.
  //
  typedef ::xsd::cxx::tree::exception< char > exception;
  typedef ::xsd::cxx::tree::bounds< char > bounds;
  typedef ::xsd::cxx::tree::duplicate_id< char > duplicate_id;
  typedef ::xsd::cxx::tree::parsing< char > parsing;
  typedef ::xsd::cxx::tree::expected_element< char > expected_element;
  typedef ::xsd::cxx::tree::unexpected_element< char > unexpected_element;
  typedef ::xsd::cxx::tree::expected_attribute< char > expected_attribute;
  typedef ::xsd::cxx::tree::unexpected_enumerator< char > unexpected_enumerator;
  typedef ::xsd::cxx::tree::expected_text_content< char > expected_text_content;
  typedef ::xsd::cxx::tree::no_prefix_mapping< char > no_prefix_mapping;

  // Error handler callback interface.
  //
  typedef ::xsd::cxx::xml::error_handler< char > error_handler;

  // DOM interaction.
  //
  namespace dom
  {
    // Automatic pointer for DOMDocument.
    //
    using ::xsd::cxx::xml::dom::auto_ptr;

#ifndef XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
#define XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
    // DOM user data key for back pointers to tree nodes.
    //
    const XMLCh* const tree_node_key = ::xsd::cxx::tree::user_data_keys::node;
#endif
  }
}

// Forward declarations.
//
namespace PSE_Molekulardynamik_WS12
{
  class nonEmptyString_t;
  class postiveDecimal_t;
  class inputType_t;
  class inputFile_t;
  class inputFiles_t;
  class potential_t;
  class simulation_t;
}


#include <memory>    // std::auto_ptr
#include <limits>    // std::numeric_limits
#include <algorithm> // std::binary_search

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/containers.hxx>
#include <xsd/cxx/tree/list.hxx>

#include <xsd/cxx/xml/dom/parsing-header.hxx>

namespace PSE_Molekulardynamik_WS12
{
  class nonEmptyString_t: public ::xml_schema::string
  {
    public:
    // Constructors.
    //
    nonEmptyString_t ();

    nonEmptyString_t (const char*);

    nonEmptyString_t (const ::std::string&);

    nonEmptyString_t (const ::xml_schema::string&);

    nonEmptyString_t (const ::xercesc::DOMElement& e,
                      ::xml_schema::flags f = 0,
                      ::xml_schema::container* c = 0);

    nonEmptyString_t (const ::xercesc::DOMAttr& a,
                      ::xml_schema::flags f = 0,
                      ::xml_schema::container* c = 0);

    nonEmptyString_t (const ::std::string& s,
                      const ::xercesc::DOMElement* e,
                      ::xml_schema::flags f = 0,
                      ::xml_schema::container* c = 0);

    nonEmptyString_t (const nonEmptyString_t& x,
                      ::xml_schema::flags f = 0,
                      ::xml_schema::container* c = 0);

    virtual nonEmptyString_t*
    _clone (::xml_schema::flags f = 0,
            ::xml_schema::container* c = 0) const;

    virtual 
    ~nonEmptyString_t ();
  };

  class postiveDecimal_t: public ::xsd::cxx::tree::fundamental_base< ::xml_schema::decimal, char, ::xml_schema::simple_type, ::xsd::cxx::tree::schema_type::decimal >
  {
    public:
    // Constructors.
    //
    postiveDecimal_t (const ::xml_schema::decimal&);

    postiveDecimal_t (const ::xercesc::DOMElement& e,
                      ::xml_schema::flags f = 0,
                      ::xml_schema::container* c = 0);

    postiveDecimal_t (const ::xercesc::DOMAttr& a,
                      ::xml_schema::flags f = 0,
                      ::xml_schema::container* c = 0);

    postiveDecimal_t (const ::std::string& s,
                      const ::xercesc::DOMElement* e,
                      ::xml_schema::flags f = 0,
                      ::xml_schema::container* c = 0);

    postiveDecimal_t (const postiveDecimal_t& x,
                      ::xml_schema::flags f = 0,
                      ::xml_schema::container* c = 0);

    virtual postiveDecimal_t*
    _clone (::xml_schema::flags f = 0,
            ::xml_schema::container* c = 0) const;

    virtual 
    ~postiveDecimal_t ();
  };

  class inputType_t: public ::xml_schema::string
  {
    public:
    enum value
    {
      list,
      cuboid
    };

    inputType_t (value v);

    inputType_t (const char* v);

    inputType_t (const ::std::string& v);

    inputType_t (const ::xml_schema::string& v);

    inputType_t (const ::xercesc::DOMElement& e,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    inputType_t (const ::xercesc::DOMAttr& a,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    inputType_t (const ::std::string& s,
                 const ::xercesc::DOMElement* e,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    inputType_t (const inputType_t& x,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    virtual inputType_t*
    _clone (::xml_schema::flags f = 0,
            ::xml_schema::container* c = 0) const;

    inputType_t&
    operator= (value v);

    virtual
    operator value () const
    {
      return _xsd_inputType_t_convert ();
    }

    protected:
    value
    _xsd_inputType_t_convert () const;

    public:
    static const char* const _xsd_inputType_t_literals_[2];
    static const value _xsd_inputType_t_indexes_[2];
  };

  class inputFile_t: public ::PSE_Molekulardynamik_WS12::nonEmptyString_t
  {
    public:
    // type
    // 
    typedef ::PSE_Molekulardynamik_WS12::inputType_t type_type;
    typedef ::xsd::cxx::tree::traits< type_type, char > type_traits;

    const type_type&
    type () const;

    type_type&
    type ();

    void
    type (const type_type& x);

    void
    type (::std::auto_ptr< type_type > p);

    // Constructors.
    //
    inputFile_t (const type_type&);

    inputFile_t (const char*,
                 const type_type&);

    inputFile_t (const ::std::string&,
                 const type_type&);

    inputFile_t (const ::xml_schema::string&,
                 const type_type&);

    inputFile_t (const ::xercesc::DOMElement& e,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    inputFile_t (const inputFile_t& x,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    virtual inputFile_t*
    _clone (::xml_schema::flags f = 0,
            ::xml_schema::container* c = 0) const;

    virtual 
    ~inputFile_t ();

    // Implementation.
    //
    protected:
    void
    parse (::xsd::cxx::xml::dom::parser< char >&,
           ::xml_schema::flags);

    protected:
    ::xsd::cxx::tree::one< type_type > type_;
  };

  class inputFiles_t: public ::xml_schema::type
  {
    public:
    // inputFile
    // 
    typedef ::PSE_Molekulardynamik_WS12::inputFile_t inputFile_type;
    typedef ::xsd::cxx::tree::sequence< inputFile_type > inputFile_sequence;
    typedef inputFile_sequence::iterator inputFile_iterator;
    typedef inputFile_sequence::const_iterator inputFile_const_iterator;
    typedef ::xsd::cxx::tree::traits< inputFile_type, char > inputFile_traits;

    const inputFile_sequence&
    inputFile () const;

    inputFile_sequence&
    inputFile ();

    void
    inputFile (const inputFile_sequence& s);

    // Constructors.
    //
    inputFiles_t ();

    inputFiles_t (const ::xercesc::DOMElement& e,
                  ::xml_schema::flags f = 0,
                  ::xml_schema::container* c = 0);

    inputFiles_t (const inputFiles_t& x,
                  ::xml_schema::flags f = 0,
                  ::xml_schema::container* c = 0);

    virtual inputFiles_t*
    _clone (::xml_schema::flags f = 0,
            ::xml_schema::container* c = 0) const;

    virtual 
    ~inputFiles_t ();

    // Implementation.
    //
    protected:
    void
    parse (::xsd::cxx::xml::dom::parser< char >&,
           ::xml_schema::flags);

    protected:
    inputFile_sequence inputFile_;
  };

  class potential_t: public ::xml_schema::string
  {
    public:
    enum value
    {
      gravitational,
      lenard_jones
    };

    potential_t (value v);

    potential_t (const char* v);

    potential_t (const ::std::string& v);

    potential_t (const ::xml_schema::string& v);

    potential_t (const ::xercesc::DOMElement& e,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    potential_t (const ::xercesc::DOMAttr& a,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    potential_t (const ::std::string& s,
                 const ::xercesc::DOMElement* e,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    potential_t (const potential_t& x,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

    virtual potential_t*
    _clone (::xml_schema::flags f = 0,
            ::xml_schema::container* c = 0) const;

    potential_t&
    operator= (value v);

    virtual
    operator value () const
    {
      return _xsd_potential_t_convert ();
    }

    protected:
    value
    _xsd_potential_t_convert () const;

    public:
    static const char* const _xsd_potential_t_literals_[2];
    static const value _xsd_potential_t_indexes_[2];
  };

  class simulation_t: public ::xml_schema::type
  {
    public:
    // outputFile
    // 
    typedef ::PSE_Molekulardynamik_WS12::nonEmptyString_t outputFile_type;
    typedef ::xsd::cxx::tree::traits< outputFile_type, char > outputFile_traits;

    const outputFile_type&
    outputFile () const;

    outputFile_type&
    outputFile ();

    void
    outputFile (const outputFile_type& x);

    void
    outputFile (::std::auto_ptr< outputFile_type > p);

    // inputFiles
    // 
    typedef ::PSE_Molekulardynamik_WS12::inputFiles_t inputFiles_type;
    typedef ::xsd::cxx::tree::traits< inputFiles_type, char > inputFiles_traits;

    const inputFiles_type&
    inputFiles () const;

    inputFiles_type&
    inputFiles ();

    void
    inputFiles (const inputFiles_type& x);

    void
    inputFiles (::std::auto_ptr< inputFiles_type > p);

    // writeFreqency
    // 
    typedef ::xml_schema::positive_integer writeFreqency_type;
    typedef ::xsd::cxx::tree::traits< writeFreqency_type, char > writeFreqency_traits;

    const writeFreqency_type&
    writeFreqency () const;

    writeFreqency_type&
    writeFreqency ();

    void
    writeFreqency (const writeFreqency_type& x);

    // t_end
    // 
    typedef ::PSE_Molekulardynamik_WS12::postiveDecimal_t t_end_type;
    typedef ::xsd::cxx::tree::traits< t_end_type, char > t_end_traits;

    const t_end_type&
    t_end () const;

    t_end_type&
    t_end ();

    void
    t_end (const t_end_type& x);

    void
    t_end (::std::auto_ptr< t_end_type > p);

    // delta_t
    // 
    typedef ::PSE_Molekulardynamik_WS12::postiveDecimal_t delta_t_type;
    typedef ::xsd::cxx::tree::traits< delta_t_type, char > delta_t_traits;

    const delta_t_type&
    delta_t () const;

    delta_t_type&
    delta_t ();

    void
    delta_t (const delta_t_type& x);

    void
    delta_t (::std::auto_ptr< delta_t_type > p);

    // potential
    // 
    typedef ::PSE_Molekulardynamik_WS12::potential_t potential_type;
    typedef ::xsd::cxx::tree::traits< potential_type, char > potential_traits;

    const potential_type&
    potential () const;

    potential_type&
    potential ();

    void
    potential (const potential_type& x);

    void
    potential (::std::auto_ptr< potential_type > p);

    // Constructors.
    //
    simulation_t (const outputFile_type&,
                  const inputFiles_type&,
                  const writeFreqency_type&,
                  const t_end_type&,
                  const delta_t_type&,
                  const potential_type&);

    simulation_t (const outputFile_type&,
                  ::std::auto_ptr< inputFiles_type >&,
                  const writeFreqency_type&,
                  const t_end_type&,
                  const delta_t_type&,
                  const potential_type&);

    simulation_t (const ::xercesc::DOMElement& e,
                  ::xml_schema::flags f = 0,
                  ::xml_schema::container* c = 0);

    simulation_t (const simulation_t& x,
                  ::xml_schema::flags f = 0,
                  ::xml_schema::container* c = 0);

    virtual simulation_t*
    _clone (::xml_schema::flags f = 0,
            ::xml_schema::container* c = 0) const;

    virtual 
    ~simulation_t ();

    // Implementation.
    //
    protected:
    void
    parse (::xsd::cxx::xml::dom::parser< char >&,
           ::xml_schema::flags);

    protected:
    ::xsd::cxx::tree::one< outputFile_type > outputFile_;
    ::xsd::cxx::tree::one< inputFiles_type > inputFiles_;
    ::xsd::cxx::tree::one< writeFreqency_type > writeFreqency_;
    ::xsd::cxx::tree::one< t_end_type > t_end_;
    ::xsd::cxx::tree::one< delta_t_type > delta_t_;
    ::xsd::cxx::tree::one< potential_type > potential_;
  };
}

#include <iosfwd>

#include <xercesc/sax/InputSource.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>

namespace PSE_Molekulardynamik_WS12
{
  // Parse a URI or a local file.
  //

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (const ::std::string& uri,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (const ::std::string& uri,
              ::xml_schema::error_handler& eh,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (const ::std::string& uri,
              ::xercesc::DOMErrorHandler& eh,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  // Parse std::istream.
  //

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              ::xml_schema::error_handler& eh,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              ::xercesc::DOMErrorHandler& eh,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              const ::std::string& id,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              const ::std::string& id,
              ::xml_schema::error_handler& eh,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              const ::std::string& id,
              ::xercesc::DOMErrorHandler& eh,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  // Parse xercesc::InputSource.
  //

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::xercesc::InputSource& is,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::xercesc::InputSource& is,
              ::xml_schema::error_handler& eh,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::xercesc::InputSource& is,
              ::xercesc::DOMErrorHandler& eh,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  // Parse xercesc::DOMDocument.
  //

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (const ::xercesc::DOMDocument& d,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument >& d,
              ::xml_schema::flags f = 0,
              const ::xml_schema::properties& p = ::xml_schema::properties ());
}

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

#endif // INPUT_PARAMETERS_HXX
