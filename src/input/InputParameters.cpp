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

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/pre.hxx>

#include "InputParameters.h"

namespace PSE_Molekulardynamik_WS12
{
  // nonEmptyString_t
  // 


  // positiveDecimal_t
  // 


  // inputType_t
  // 

  inputType_t::
  inputType_t (value v)
  : ::xml_schema::string (_xsd_inputType_t_literals_[v])
  {
  }

  inputType_t::
  inputType_t (const char* v)
  : ::xml_schema::string (v)
  {
  }

  inputType_t::
  inputType_t (const ::std::string& v)
  : ::xml_schema::string (v)
  {
  }

  inputType_t::
  inputType_t (const ::xml_schema::string& v)
  : ::xml_schema::string (v)
  {
  }

  inputType_t::
  inputType_t (const inputType_t& v,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::xml_schema::string (v, f, c)
  {
  }

  inputType_t& inputType_t::
  operator= (value v)
  {
    static_cast< ::xml_schema::string& > (*this) = 
    ::xml_schema::string (_xsd_inputType_t_literals_[v]);

    return *this;
  }


  // inputFile_t
  // 

  const inputFile_t::type_type& inputFile_t::
  type () const
  {
    return this->type_.get ();
  }

  inputFile_t::type_type& inputFile_t::
  type ()
  {
    return this->type_.get ();
  }

  void inputFile_t::
  type (const type_type& x)
  {
    this->type_.set (x);
  }

  void inputFile_t::
  type (::std::auto_ptr< type_type > x)
  {
    this->type_.set (x);
  }


  // inputFiles_t
  // 

  const inputFiles_t::inputFile_sequence& inputFiles_t::
  inputFile () const
  {
    return this->inputFile_;
  }

  inputFiles_t::inputFile_sequence& inputFiles_t::
  inputFile ()
  {
    return this->inputFile_;
  }

  void inputFiles_t::
  inputFile (const inputFile_sequence& s)
  {
    this->inputFile_ = s;
  }


  // potential_t
  // 

  potential_t::
  potential_t (value v)
  : ::xml_schema::string (_xsd_potential_t_literals_[v])
  {
  }

  potential_t::
  potential_t (const char* v)
  : ::xml_schema::string (v)
  {
  }

  potential_t::
  potential_t (const ::std::string& v)
  : ::xml_schema::string (v)
  {
  }

  potential_t::
  potential_t (const ::xml_schema::string& v)
  : ::xml_schema::string (v)
  {
  }

  potential_t::
  potential_t (const potential_t& v,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::xml_schema::string (v, f, c)
  {
  }

  potential_t& potential_t::
  operator= (value v)
  {
    static_cast< ::xml_schema::string& > (*this) = 
    ::xml_schema::string (_xsd_potential_t_literals_[v]);

    return *this;
  }


  // simulation_t
  // 

  const simulation_t::outputFile_type& simulation_t::
  outputFile () const
  {
    return this->outputFile_.get ();
  }

  simulation_t::outputFile_type& simulation_t::
  outputFile ()
  {
    return this->outputFile_.get ();
  }

  void simulation_t::
  outputFile (const outputFile_type& x)
  {
    this->outputFile_.set (x);
  }

  void simulation_t::
  outputFile (::std::auto_ptr< outputFile_type > x)
  {
    this->outputFile_.set (x);
  }

  const simulation_t::inputFiles_type& simulation_t::
  inputFiles () const
  {
    return this->inputFiles_.get ();
  }

  simulation_t::inputFiles_type& simulation_t::
  inputFiles ()
  {
    return this->inputFiles_.get ();
  }

  void simulation_t::
  inputFiles (const inputFiles_type& x)
  {
    this->inputFiles_.set (x);
  }

  void simulation_t::
  inputFiles (::std::auto_ptr< inputFiles_type > x)
  {
    this->inputFiles_.set (x);
  }

  const simulation_t::writeFrequency_type& simulation_t::
  writeFrequency () const
  {
    return this->writeFrequency_.get ();
  }

  simulation_t::writeFrequency_type& simulation_t::
  writeFrequency ()
  {
    return this->writeFrequency_.get ();
  }

  void simulation_t::
  writeFrequency (const writeFrequency_type& x)
  {
    this->writeFrequency_.set (x);
  }

  const simulation_t::t_end_type& simulation_t::
  t_end () const
  {
    return this->t_end_.get ();
  }

  simulation_t::t_end_type& simulation_t::
  t_end ()
  {
    return this->t_end_.get ();
  }

  void simulation_t::
  t_end (const t_end_type& x)
  {
    this->t_end_.set (x);
  }

  void simulation_t::
  t_end (::std::auto_ptr< t_end_type > x)
  {
    this->t_end_.set (x);
  }

  const simulation_t::delta_t_type& simulation_t::
  delta_t () const
  {
    return this->delta_t_.get ();
  }

  simulation_t::delta_t_type& simulation_t::
  delta_t ()
  {
    return this->delta_t_.get ();
  }

  void simulation_t::
  delta_t (const delta_t_type& x)
  {
    this->delta_t_.set (x);
  }

  void simulation_t::
  delta_t (::std::auto_ptr< delta_t_type > x)
  {
    this->delta_t_.set (x);
  }

  const simulation_t::potential_type& simulation_t::
  potential () const
  {
    return this->potential_.get ();
  }

  simulation_t::potential_type& simulation_t::
  potential ()
  {
    return this->potential_.get ();
  }

  void simulation_t::
  potential (const potential_type& x)
  {
    this->potential_.set (x);
  }

  void simulation_t::
  potential (::std::auto_ptr< potential_type > x)
  {
    this->potential_.set (x);
  }
}

#include <xsd/cxx/xml/dom/parsing-source.hxx>

namespace PSE_Molekulardynamik_WS12
{
  // nonEmptyString_t
  //

  nonEmptyString_t::
  nonEmptyString_t ()
  : ::xml_schema::string ()
  {
  }

  nonEmptyString_t::
  nonEmptyString_t (const char* _xsd_string_base)
  : ::xml_schema::string (_xsd_string_base)
  {
  }

  nonEmptyString_t::
  nonEmptyString_t (const ::std::string& _xsd_string_base)
  : ::xml_schema::string (_xsd_string_base)
  {
  }

  nonEmptyString_t::
  nonEmptyString_t (const ::xml_schema::string& _xsd_string_base)
  : ::xml_schema::string (_xsd_string_base)
  {
  }

  nonEmptyString_t::
  nonEmptyString_t (const nonEmptyString_t& x,
                    ::xml_schema::flags f,
                    ::xml_schema::container* c)
  : ::xml_schema::string (x, f, c)
  {
  }

  nonEmptyString_t::
  nonEmptyString_t (const ::xercesc::DOMElement& e,
                    ::xml_schema::flags f,
                    ::xml_schema::container* c)
  : ::xml_schema::string (e, f, c)
  {
  }

  nonEmptyString_t::
  nonEmptyString_t (const ::xercesc::DOMAttr& a,
                    ::xml_schema::flags f,
                    ::xml_schema::container* c)
  : ::xml_schema::string (a, f, c)
  {
  }

  nonEmptyString_t::
  nonEmptyString_t (const ::std::string& s,
                    const ::xercesc::DOMElement* e,
                    ::xml_schema::flags f,
                    ::xml_schema::container* c)
  : ::xml_schema::string (s, e, f, c)
  {
  }

  nonEmptyString_t* nonEmptyString_t::
  _clone (::xml_schema::flags f,
          ::xml_schema::container* c) const
  {
    return new class nonEmptyString_t (*this, f, c);
  }

  nonEmptyString_t::
  ~nonEmptyString_t ()
  {
  }

  // positiveDecimal_t
  //

  positiveDecimal_t::
  positiveDecimal_t (const ::xml_schema::decimal& _xsd_decimal_base)
  : ::xsd::cxx::tree::fundamental_base< ::xml_schema::decimal, char, ::xml_schema::simple_type, ::xsd::cxx::tree::schema_type::decimal > (_xsd_decimal_base)
  {
  }

  positiveDecimal_t::
  positiveDecimal_t (const positiveDecimal_t& x,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
  : ::xsd::cxx::tree::fundamental_base< ::xml_schema::decimal, char, ::xml_schema::simple_type, ::xsd::cxx::tree::schema_type::decimal > (x, f, c)
  {
  }

  positiveDecimal_t::
  positiveDecimal_t (const ::xercesc::DOMElement& e,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
  : ::xsd::cxx::tree::fundamental_base< ::xml_schema::decimal, char, ::xml_schema::simple_type, ::xsd::cxx::tree::schema_type::decimal > (e, f, c)
  {
  }

  positiveDecimal_t::
  positiveDecimal_t (const ::xercesc::DOMAttr& a,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
  : ::xsd::cxx::tree::fundamental_base< ::xml_schema::decimal, char, ::xml_schema::simple_type, ::xsd::cxx::tree::schema_type::decimal > (a, f, c)
  {
  }

  positiveDecimal_t::
  positiveDecimal_t (const ::std::string& s,
                     const ::xercesc::DOMElement* e,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
  : ::xsd::cxx::tree::fundamental_base< ::xml_schema::decimal, char, ::xml_schema::simple_type, ::xsd::cxx::tree::schema_type::decimal > (s, e, f, c)
  {
  }

  positiveDecimal_t* positiveDecimal_t::
  _clone (::xml_schema::flags f,
          ::xml_schema::container* c) const
  {
    return new class positiveDecimal_t (*this, f, c);
  }

  positiveDecimal_t::
  ~positiveDecimal_t ()
  {
  }

  // inputType_t
  //

  inputType_t::
  inputType_t (const ::xercesc::DOMElement& e,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::xml_schema::string (e, f, c)
  {
    _xsd_inputType_t_convert ();
  }

  inputType_t::
  inputType_t (const ::xercesc::DOMAttr& a,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::xml_schema::string (a, f, c)
  {
    _xsd_inputType_t_convert ();
  }

  inputType_t::
  inputType_t (const ::std::string& s,
               const ::xercesc::DOMElement* e,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::xml_schema::string (s, e, f, c)
  {
    _xsd_inputType_t_convert ();
  }

  inputType_t* inputType_t::
  _clone (::xml_schema::flags f,
          ::xml_schema::container* c) const
  {
    return new class inputType_t (*this, f, c);
  }

  inputType_t::value inputType_t::
  _xsd_inputType_t_convert () const
  {
    ::xsd::cxx::tree::enum_comparator< char > c (_xsd_inputType_t_literals_);
    const value* i (::std::lower_bound (
                      _xsd_inputType_t_indexes_,
                      _xsd_inputType_t_indexes_ + 2,
                      *this,
                      c));

    if (i == _xsd_inputType_t_indexes_ + 2 || _xsd_inputType_t_literals_[*i] != *this)
    {
      throw ::xsd::cxx::tree::unexpected_enumerator < char > (*this);
    }

    return *i;
  }

  const char* const inputType_t::
  _xsd_inputType_t_literals_[2] =
  {
    "list",
    "cuboid"
  };

  const inputType_t::value inputType_t::
  _xsd_inputType_t_indexes_[2] =
  {
    ::PSE_Molekulardynamik_WS12::inputType_t::cuboid,
    ::PSE_Molekulardynamik_WS12::inputType_t::list
  };

  // inputFile_t
  //

  inputFile_t::
  inputFile_t (const type_type& type)
  : ::PSE_Molekulardynamik_WS12::nonEmptyString_t (),
    type_ (type, ::xml_schema::flags (), this)
  {
  }

  inputFile_t::
  inputFile_t (const char* _xsd_string_base,
               const type_type& type)
  : ::PSE_Molekulardynamik_WS12::nonEmptyString_t (_xsd_string_base),
    type_ (type, ::xml_schema::flags (), this)
  {
  }

  inputFile_t::
  inputFile_t (const ::std::string& _xsd_string_base,
               const type_type& type)
  : ::PSE_Molekulardynamik_WS12::nonEmptyString_t (_xsd_string_base),
    type_ (type, ::xml_schema::flags (), this)
  {
  }

  inputFile_t::
  inputFile_t (const ::xml_schema::string& _xsd_string_base,
               const type_type& type)
  : ::PSE_Molekulardynamik_WS12::nonEmptyString_t (_xsd_string_base),
    type_ (type, ::xml_schema::flags (), this)
  {
  }

  inputFile_t::
  inputFile_t (const inputFile_t& x,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::PSE_Molekulardynamik_WS12::nonEmptyString_t (x, f, c),
    type_ (x.type_, f, this)
  {
  }

  inputFile_t::
  inputFile_t (const ::xercesc::DOMElement& e,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::PSE_Molekulardynamik_WS12::nonEmptyString_t (e, f | ::xml_schema::flags::base, c),
    type_ (f, this)
  {
    if ((f & ::xml_schema::flags::base) == 0)
    {
      ::xsd::cxx::xml::dom::parser< char > p (e, false, true);
      this->parse (p, f);
    }
  }

  void inputFile_t::
  parse (::xsd::cxx::xml::dom::parser< char >& p,
         ::xml_schema::flags f)
  {
    while (p.more_attributes ())
    {
      const ::xercesc::DOMAttr& i (p.next_attribute ());
      const ::xsd::cxx::xml::qualified_name< char > n (
        ::xsd::cxx::xml::dom::name< char > (i));

      if (n.name () == "type" && n.namespace_ ().empty ())
      {
        ::std::auto_ptr< type_type > r (
          type_traits::create (i, f, this));

        this->type_.set (r);
        continue;
      }
    }

    if (!type_.present ())
    {
      throw ::xsd::cxx::tree::expected_attribute< char > (
        "type",
        "");
    }
  }

  inputFile_t* inputFile_t::
  _clone (::xml_schema::flags f,
          ::xml_schema::container* c) const
  {
    return new class inputFile_t (*this, f, c);
  }

  inputFile_t::
  ~inputFile_t ()
  {
  }

  // inputFiles_t
  //

  inputFiles_t::
  inputFiles_t ()
  : ::xml_schema::type (),
    inputFile_ (::xml_schema::flags (), this)
  {
  }

  inputFiles_t::
  inputFiles_t (const inputFiles_t& x,
                ::xml_schema::flags f,
                ::xml_schema::container* c)
  : ::xml_schema::type (x, f, c),
    inputFile_ (x.inputFile_, f, this)
  {
  }

  inputFiles_t::
  inputFiles_t (const ::xercesc::DOMElement& e,
                ::xml_schema::flags f,
                ::xml_schema::container* c)
  : ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
    inputFile_ (f, this)
  {
    if ((f & ::xml_schema::flags::base) == 0)
    {
      ::xsd::cxx::xml::dom::parser< char > p (e, true, false);
      this->parse (p, f);
    }
  }

  void inputFiles_t::
  parse (::xsd::cxx::xml::dom::parser< char >& p,
         ::xml_schema::flags f)
  {
    for (; p.more_elements (); p.next_element ())
    {
      const ::xercesc::DOMElement& i (p.cur_element ());
      const ::xsd::cxx::xml::qualified_name< char > n (
        ::xsd::cxx::xml::dom::name< char > (i));

      // inputFile
      //
      if (n.name () == "inputFile" && n.namespace_ () == "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12")
      {
        ::std::auto_ptr< inputFile_type > r (
          inputFile_traits::create (i, f, this));

        this->inputFile_.push_back (r);
        continue;
      }

      break;
    }
  }

  inputFiles_t* inputFiles_t::
  _clone (::xml_schema::flags f,
          ::xml_schema::container* c) const
  {
    return new class inputFiles_t (*this, f, c);
  }

  inputFiles_t::
  ~inputFiles_t ()
  {
  }

  // potential_t
  //

  potential_t::
  potential_t (const ::xercesc::DOMElement& e,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::xml_schema::string (e, f, c)
  {
    _xsd_potential_t_convert ();
  }

  potential_t::
  potential_t (const ::xercesc::DOMAttr& a,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::xml_schema::string (a, f, c)
  {
    _xsd_potential_t_convert ();
  }

  potential_t::
  potential_t (const ::std::string& s,
               const ::xercesc::DOMElement* e,
               ::xml_schema::flags f,
               ::xml_schema::container* c)
  : ::xml_schema::string (s, e, f, c)
  {
    _xsd_potential_t_convert ();
  }

  potential_t* potential_t::
  _clone (::xml_schema::flags f,
          ::xml_schema::container* c) const
  {
    return new class potential_t (*this, f, c);
  }

  potential_t::value potential_t::
  _xsd_potential_t_convert () const
  {
    ::xsd::cxx::tree::enum_comparator< char > c (_xsd_potential_t_literals_);
    const value* i (::std::lower_bound (
                      _xsd_potential_t_indexes_,
                      _xsd_potential_t_indexes_ + 2,
                      *this,
                      c));

    if (i == _xsd_potential_t_indexes_ + 2 || _xsd_potential_t_literals_[*i] != *this)
    {
      throw ::xsd::cxx::tree::unexpected_enumerator < char > (*this);
    }

    return *i;
  }

  const char* const potential_t::
  _xsd_potential_t_literals_[2] =
  {
    "gravitational",
    "lenard jones"
  };

  const potential_t::value potential_t::
  _xsd_potential_t_indexes_[2] =
  {
    ::PSE_Molekulardynamik_WS12::potential_t::gravitational,
    ::PSE_Molekulardynamik_WS12::potential_t::lenard_jones
  };

  // simulation_t
  //

  simulation_t::
  simulation_t (const outputFile_type& outputFile,
                const inputFiles_type& inputFiles,
                const writeFrequency_type& writeFrequency,
                const t_end_type& t_end,
                const delta_t_type& delta_t,
                const potential_type& potential)
  : ::xml_schema::type (),
    outputFile_ (outputFile, ::xml_schema::flags (), this),
    inputFiles_ (inputFiles, ::xml_schema::flags (), this),
    writeFrequency_ (writeFrequency, ::xml_schema::flags (), this),
    t_end_ (t_end, ::xml_schema::flags (), this),
    delta_t_ (delta_t, ::xml_schema::flags (), this),
    potential_ (potential, ::xml_schema::flags (), this)
  {
  }

  simulation_t::
  simulation_t (const outputFile_type& outputFile,
                ::std::auto_ptr< inputFiles_type >& inputFiles,
                const writeFrequency_type& writeFrequency,
                const t_end_type& t_end,
                const delta_t_type& delta_t,
                const potential_type& potential)
  : ::xml_schema::type (),
    outputFile_ (outputFile, ::xml_schema::flags (), this),
    inputFiles_ (inputFiles, ::xml_schema::flags (), this),
    writeFrequency_ (writeFrequency, ::xml_schema::flags (), this),
    t_end_ (t_end, ::xml_schema::flags (), this),
    delta_t_ (delta_t, ::xml_schema::flags (), this),
    potential_ (potential, ::xml_schema::flags (), this)
  {
  }

  simulation_t::
  simulation_t (const simulation_t& x,
                ::xml_schema::flags f,
                ::xml_schema::container* c)
  : ::xml_schema::type (x, f, c),
    outputFile_ (x.outputFile_, f, this),
    inputFiles_ (x.inputFiles_, f, this),
    writeFrequency_ (x.writeFrequency_, f, this),
    t_end_ (x.t_end_, f, this),
    delta_t_ (x.delta_t_, f, this),
    potential_ (x.potential_, f, this)
  {
  }

  simulation_t::
  simulation_t (const ::xercesc::DOMElement& e,
                ::xml_schema::flags f,
                ::xml_schema::container* c)
  : ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
    outputFile_ (f, this),
    inputFiles_ (f, this),
    writeFrequency_ (f, this),
    t_end_ (f, this),
    delta_t_ (f, this),
    potential_ (f, this)
  {
    if ((f & ::xml_schema::flags::base) == 0)
    {
      ::xsd::cxx::xml::dom::parser< char > p (e, true, false);
      this->parse (p, f);
    }
  }

  void simulation_t::
  parse (::xsd::cxx::xml::dom::parser< char >& p,
         ::xml_schema::flags f)
  {
    for (; p.more_elements (); p.next_element ())
    {
      const ::xercesc::DOMElement& i (p.cur_element ());
      const ::xsd::cxx::xml::qualified_name< char > n (
        ::xsd::cxx::xml::dom::name< char > (i));

      // outputFile
      //
      if (n.name () == "outputFile" && n.namespace_ () == "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12")
      {
        ::std::auto_ptr< outputFile_type > r (
          outputFile_traits::create (i, f, this));

        if (!outputFile_.present ())
        {
          this->outputFile_.set (r);
          continue;
        }
      }

      // inputFiles
      //
      if (n.name () == "inputFiles" && n.namespace_ () == "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12")
      {
        ::std::auto_ptr< inputFiles_type > r (
          inputFiles_traits::create (i, f, this));

        if (!inputFiles_.present ())
        {
          this->inputFiles_.set (r);
          continue;
        }
      }

      // writeFrequency
      //
      if (n.name () == "writeFrequency" && n.namespace_ () == "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12")
      {
        if (!writeFrequency_.present ())
        {
          this->writeFrequency_.set (writeFrequency_traits::create (i, f, this));
          continue;
        }
      }

      // t_end
      //
      if (n.name () == "t_end" && n.namespace_ () == "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12")
      {
        ::std::auto_ptr< t_end_type > r (
          t_end_traits::create (i, f, this));

        if (!t_end_.present ())
        {
          this->t_end_.set (r);
          continue;
        }
      }

      // delta_t
      //
      if (n.name () == "delta_t" && n.namespace_ () == "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12")
      {
        ::std::auto_ptr< delta_t_type > r (
          delta_t_traits::create (i, f, this));

        if (!delta_t_.present ())
        {
          this->delta_t_.set (r);
          continue;
        }
      }

      // potential
      //
      if (n.name () == "potential" && n.namespace_ () == "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12")
      {
        ::std::auto_ptr< potential_type > r (
          potential_traits::create (i, f, this));

        if (!potential_.present ())
        {
          this->potential_.set (r);
          continue;
        }
      }

      break;
    }

    if (!outputFile_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "outputFile",
        "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12");
    }

    if (!inputFiles_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "inputFiles",
        "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12");
    }

    if (!writeFrequency_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "writeFrequency",
        "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12");
    }

    if (!t_end_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "t_end",
        "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12");
    }

    if (!delta_t_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "delta_t",
        "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12");
    }

    if (!potential_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "potential",
        "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12");
    }
  }

  simulation_t* simulation_t::
  _clone (::xml_schema::flags f,
          ::xml_schema::container* c) const
  {
    return new class simulation_t (*this, f, c);
  }

  simulation_t::
  ~simulation_t ()
  {
  }
}

#include <istream>
#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/tree/error-handler.hxx>

namespace PSE_Molekulardynamik_WS12
{
  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (const ::std::string& u,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::tree::error_handler< char > h;

    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        u, h, p, f));

    h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

    ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t > r (
      ::PSE_Molekulardynamik_WS12::simulation (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (const ::std::string& u,
              ::xml_schema::error_handler& h,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        u, h, p, f));

    if (!d.get ())
      throw ::xsd::cxx::tree::parsing< char > ();

    ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t > r (
      ::PSE_Molekulardynamik_WS12::simulation (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (const ::std::string& u,
              ::xercesc::DOMErrorHandler& h,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        u, h, p, f));

    if (!d.get ())
      throw ::xsd::cxx::tree::parsing< char > ();

    ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t > r (
      ::PSE_Molekulardynamik_WS12::simulation (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::xml::sax::std_input_source isrc (is);
    return ::PSE_Molekulardynamik_WS12::simulation (isrc, f, p);
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              ::xml_schema::error_handler& h,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::xml::sax::std_input_source isrc (is);
    return ::PSE_Molekulardynamik_WS12::simulation (isrc, h, f, p);
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              ::xercesc::DOMErrorHandler& h,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::sax::std_input_source isrc (is);
    return ::PSE_Molekulardynamik_WS12::simulation (isrc, h, f, p);
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              const ::std::string& sid,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
    return ::PSE_Molekulardynamik_WS12::simulation (isrc, f, p);
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              const ::std::string& sid,
              ::xml_schema::error_handler& h,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
    return ::PSE_Molekulardynamik_WS12::simulation (isrc, h, f, p);
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::std::istream& is,
              const ::std::string& sid,
              ::xercesc::DOMErrorHandler& h,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
    return ::PSE_Molekulardynamik_WS12::simulation (isrc, h, f, p);
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::xercesc::InputSource& i,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xsd::cxx::tree::error_handler< char > h;

    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        i, h, p, f));

    h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

    ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t > r (
      ::PSE_Molekulardynamik_WS12::simulation (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::xercesc::InputSource& i,
              ::xml_schema::error_handler& h,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        i, h, p, f));

    if (!d.get ())
      throw ::xsd::cxx::tree::parsing< char > ();

    ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t > r (
      ::PSE_Molekulardynamik_WS12::simulation (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::xercesc::InputSource& i,
              ::xercesc::DOMErrorHandler& h,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        i, h, p, f));

    if (!d.get ())
      throw ::xsd::cxx::tree::parsing< char > ();

    ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t > r (
      ::PSE_Molekulardynamik_WS12::simulation (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (const ::xercesc::DOMDocument& d,
              ::xml_schema::flags f,
              const ::xml_schema::properties& p)
  {
    if (f & ::xml_schema::flags::keep_dom)
    {
      ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > c (
        static_cast< ::xercesc::DOMDocument* > (d.cloneNode (true)));

      ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t > r (
        ::PSE_Molekulardynamik_WS12::simulation (
          c, f | ::xml_schema::flags::own_dom, p));

      return r;
    }

    const ::xercesc::DOMElement& e (*d.getDocumentElement ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (e));

    if (n.name () == "simulation" &&
        n.namespace_ () == "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12")
    {
      ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t > r (
        ::xsd::cxx::tree::traits< ::PSE_Molekulardynamik_WS12::simulation_t, char >::create (
          e, f, 0));
      return r;
    }

    throw ::xsd::cxx::tree::unexpected_element < char > (
      n.name (),
      n.namespace_ (),
      "simulation",
      "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12");
  }

  ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t >
  simulation (::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument >& d,
              ::xml_schema::flags f,
              const ::xml_schema::properties&)
  {
    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > c (
      ((f & ::xml_schema::flags::keep_dom) &&
       !(f & ::xml_schema::flags::own_dom))
      ? static_cast< ::xercesc::DOMDocument* > (d->cloneNode (true))
      : 0);

    ::xercesc::DOMDocument& doc (c.get () ? *c : *d);
    const ::xercesc::DOMElement& e (*doc.getDocumentElement ());

    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (e));

    if (f & ::xml_schema::flags::keep_dom)
      doc.setUserData (::xml_schema::dom::tree_node_key,
                       (c.get () ? &c : &d),
                       0);

    if (n.name () == "simulation" &&
        n.namespace_ () == "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12")
    {
      ::std::auto_ptr< ::PSE_Molekulardynamik_WS12::simulation_t > r (
        ::xsd::cxx::tree::traits< ::PSE_Molekulardynamik_WS12::simulation_t, char >::create (
          e, f, 0));
      return r;
    }

    throw ::xsd::cxx::tree::unexpected_element < char > (
      n.name (),
      n.namespace_ (),
      "simulation",
      "http://www5.in.tum.de/wiki/index.php/PSE_Molekulardynamik_WS12");
  }
}

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

