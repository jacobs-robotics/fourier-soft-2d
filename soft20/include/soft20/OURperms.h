/***************************************************************************
  **************************************************************************
    
  Spherical Harmonic Transform Kit 2.7
    
  Copyright 1997-2003  Sean Moore, Dennis Healy,
                       Dan Rockmore, Peter Kostelec
  Copyright 2004  Peter Kostelec, Dan Rockmore

  This file is part of SpharmonicKit.

  SpharmonicKit is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  SpharmonicKit is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  See the accompanying LICENSE file for details.

  ************************************************************************
  ************************************************************************/

/*
  function to interface with the OUR permutation arrays
  for FCTs defined in OURperms.c
  */

#ifndef _OURPERMS_H
#define _OURPERMS_H

#ifdef __cplusplus
extern "C"
{
#endif

extern const int *get_perm( int ) ;

#ifdef __cplusplus
}
#endif

#endif /* _OURPERMS_H */
