/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/


#include "itkVectorImage.h"
#include "itkSLICImageFilter.h"

int itkSLICImageFilterTest(int argc, char *arg[])
{
  const unsigned int Dimension = 3;
  typedef itk::VectorImage<float, Dimension> InputImageType;
  typedef itk::Image<unsigned int, Dimension> OutputImageType;

  typedef itk::SLICImageFilter< InputImageType, OutputImageType > FilterType;

  FilterType::Pointer filter = FilterType::New();

  return 0;
}
