/*
    GPT - Grid Python Toolkit
    Copyright (C) 2020  Christoph Lehner (christoph.lehner@ur.de, https://github.com/lehner/gpt)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    Description:  We need to fail gracefully since we also run in an interpreter; infrastructure goes here
*/
#define STRX(x) #x
#define STR(x) STRX(x)
#define ASSERT(x)				\
  { if ( !(x) )throw "Assert " #x " failed in file " __FILE__ ":"  STR(__LINE__); };
#define ERR(...)							\
  { char msg[1024]; snprintf(msg,sizeof(msg)-100,__VA_ARGS__);		\
    strcat(msg, " in file " __FILE__ ":"  STR(__LINE__)); throw (const char*)msg; };

#define EXPORT(name,...)					   \
  PyObject* cgpt_ ## name(PyObject* self, PyObject* args) {	   \
    try {							   \
      __VA_ARGS__;						   \
      return NULL;						   \
    } catch (const char* err) {					   \
      PyErr_SetString(PyExc_RuntimeError,err);			   \
      return NULL;						   \
    }								   \
  }
