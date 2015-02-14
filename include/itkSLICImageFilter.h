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
#ifndef itkSLICImageFilter_h
#define itkSLICImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVariableLengthVector.h"

#include "itkShrinkImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

template< typename TInputImage, typename TOutputImage >
class SLICImageFilter:
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef SLICImageFilter                                 Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ExtractImageFilter, ImageToImageFilter);

  /** Image type information. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  typedef double DistanceImagePixelType;
  typedef Image<DistanceImagePixelType, ImageDimension> DistanceImageType;

  typedef FixedArray< unsigned int, ImageDimension > SuperGridSizeType;

  itkSetMacro(SuperGridSize, SuperGridSizeType);
  void SetSuperGridSize(unsigned int factor);
  void SetSuperGridSize(unsigned int i, unsigned int factor);


protected:
  SLICImageFilter()
    {
      m_SuperGridSize.Fill(10);
    }
  ~SLICImageFilter() {}

  void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE
    {
      Superclass::PrintSelf(os, indent);
    }


  /** Generate full output and require full input */
  void EnlargeOutputRequestedRegion(DataObject *output)
    {
      Superclass::EnlargeOutputRequestedRegion(output);
      output->SetRequestedRegionToLargestPossibleRegion();
    }



  void GenerateData() ITK_OVERRIDE
    {
    }

private:
  SLICImageFilter(const Self &);    //purposely not implemented
  void operator=(const Self &);     //purposely not implemented

  SuperGridSizeType m_SuperGridSize;

};
} // end namespace itk

#endif //itkSLICImageFilter_h
