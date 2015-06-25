/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: VTKMappedIdFilterT.h,v $
  Language:  C++
  Date:      $Date: 2002/02/01 18:13:36 $
  Version:   $Revision: 1.1 $


Copyright (c) 1993-2001 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
// .NAME vtkMappedIdFilterT - generate scalars or field data from point and cell ids
// .SECTION Description
// vtkIdFilter is a filter to that generates scalars or field data
// using cell and point ids. That is, the point attribute data scalars
// or field data are generated from the point ids, and the cell
// attribute data scalars or field data are generated from the the
// cell ids.
//
// Typically this filter is used with vtkLabeledDataMapper (and possibly
// vtkSelectVisiblePoints) to create labels for points and cells, or labels
// for the point or cell data scalar values.
//
// This filter applies a node/cell map to generate the id's
//

#ifndef __vtkMappedIdFilterT_h
#define __vtkMappedIdFilterT_h

#include "vtkIdFilter.h"
#include "vtkDataSetToDataSetFilter.h"

class VTK_GRAPHICS_EXPORT vtkMappedIdFilterT : public vtkIdFilter
{
public:
  vtkTypeMacro(vtkMappedIdFilterT,vtkIdFilter);
  //void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct object with PointIds and CellIds on; and ids being generated
  // as scalars.
  static vtkMappedIdFilterT *New();

	// set/get the point numbering map
	void SetPointMap(int* map) { PointNumberMap = map; };
	int* GetPointMap(void) const { return PointNumberMap; };

	// set/get the point numbering map
	void SetCellMap(int* map) { CellNumberMap = map; };
	int* GetCellMap(void) const { return CellNumberMap; };

protected:
  vtkMappedIdFilterT();
//  ~vtkMappedIdFilterT();

  void Execute();
  
  //point numbering map
  int* PointNumberMap;
  
  //cell numbering map
  int* CellNumberMap;

private:
  vtkMappedIdFilterT(const vtkMappedIdFilterT&);  // Not implemented.
  void operator=(const vtkMappedIdFilterT&);  // Not implemented.
};

#endif


